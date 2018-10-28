//! Numpress utility.
//!
//! A pure rust implementation of [`ms-numpress`], a fast,
//! minimally lossy compression algorithm for mass spectrometry data.
//!
//! # Additional Information
//!
//! The API makes extensive use of unsafe Rust features, and therefore
//! cannot guarantee low-level memory safety. Use at your own risk.
//!
//! [`ms-numpress`]: https://github.com/ms-numpress/ms-numpress

use std::error::Error as StdError;
use std::fmt;
use std::result::Result as StdResult;

#[cfg(test)]
extern crate rand;

#[cfg(test)]
#[macro_use]
extern crate assert_approx_eq;

// ERROR

/// Type of error encountered during compression or decompression.
#[derive(Debug)]
pub enum ErrorKind {
    /// Encoding or decoding fails due to corrupt input data.
    CorruptInputData,
    /// Next number to compress overflows `i64` type.
    OverflowError,
    /// Number is out-of-range of `[i32::min_value(), i32::max_value()]`.
    OutOfRange,
}

/// Custom error for Numpress compression.
#[derive(Debug)]
pub struct Error(ErrorKind);

impl From<ErrorKind> for Error {
    fn from(kind: ErrorKind) -> Self {
        Error(kind)
    }
}

impl Error {
    /// Get error type.
    pub fn kind(&self) -> &ErrorKind {
        &self.0
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Numpress error: {}", self.description())
    }
}

impl StdError for Error {
    fn description(&self) -> &str {
        match self.kind() {
            ErrorKind::CorruptInputData => "corrupt input data.",
            ErrorKind::OverflowError => "next number overflows.",
            ErrorKind::OutOfRange => "cannot encode number outside of [i32::min_value(), i32::max_value()].",
        }
    }

    fn cause(&self) -> Option<&StdError> {
        None
    }
}

/// Specialized result for Numpress operations.
pub type Result<T> = StdResult<T, Error>;

pub mod low_level {

use std::mem::{uninitialized, transmute};
use super::{ErrorKind, Result};

// FIXED POINT

/// Encode fixed point, disabling all spurious bounds checking for performance.
///
/// Encodes the fixed point into `dst`.
///
/// * `fixed_point`     - Floating point to encode.
/// * `dst`             - Destination buffer with at least 8 elements.
pub(crate) unsafe extern "C" fn encode_fixed_point(
    fixed_point: f64,
    dst: *mut u8
)
{
    let fp: [u8; 8] = transmute(fixed_point);
    for i in 0..8 {
        #[cfg(target_endian = "big")] {
            *dst.add(i) = *fp.get_unchecked(i);
        }
        #[cfg(target_endian = "little")] {
            *dst.add(i) = *fp.get_unchecked(7-i);
        }
    }
}

/// Decode fixed point, disabling all spurious bounds checking for performance.
///
/// Decodes the data into an f64.
///
/// * `data`            - Input buffer with at least 8 elements.
pub(crate) unsafe extern "C" fn decode_fixed_point(
    data: *const u8
)
    -> f64
{
    let mut fp: [u8; 8] = uninitialized();
    for i in 0..8 {
        #[cfg(target_endian = "big")] {
            *fp.get_unchecked_mut(i) = *data.add(i);
        }
        #[cfg(target_endian = "little")] {
            *fp.get_unchecked_mut(i) = *data.add(7-i);
        }
    }

    transmute(fp)
}

// INT

/// Encodes the int x as a number of half bytes in res.
/// res_length is incremented by the number of half bytes,
/// which will be 1 <= n <= 9.
///
/// This code cannot overflow, due to the use of multiplication in values
/// <= 8, only right bit shifts (>>).
pub(crate) unsafe extern "C" fn encode_int(
    x: u32,
    res: *mut u8,
    res_length: *mut usize
)
{
    // get the bit pattern of input
    let mut m: u32;

    const MASK: u32 = 0xf0000000;
    let init = x & MASK;

    if init == 0 {
        let mut l: u8 = 8;
        for i in 0..8 {
            m = MASK >> (4*i);
            if (x & m) != 0 {
                l = i;
                break;
            }
        }
        *res = l;
        for i in l..8 {
            let xi = x >> (4*(i-l));
            let ri = res.add(1 + (i as usize) - (l as usize));
            *ri = xi as u8;
        }
        *res_length += 1 + 8 - (l as usize);
    } else if init == MASK {
        let mut l: u8 = 7;
        for i in 0..8 {
            m = MASK >> (4*i);
            if (x & m) != m {
                l = i;
                break;
            }
        }
        *res = l + 8;
        for i in l..8 {
            let xi = x >> (4*(i-l));
            let ri = res.add(1 + (i as usize) - (l as usize));
            *ri = xi as u8;
        }
        *res_length += 1 + 8 - (l as usize);

    } else {
        *res = 0;
        for i in 0..8 {
            let xi = x >> (4*i);
            *res.add(i+1) = xi as u8;
        }
        *res_length += 9;
    }
}

/// Decodes an int from the half bytes in bp.
///
/// Lossless reverse of `encode_int`.
pub(crate) unsafe extern "C" fn decode_int(
    data: *const u8,
    di: *mut usize,
    max_di: usize,
    half: *mut usize,
    res: *mut u32
)
    -> Result<()>
{
    let n: usize;
    let mask: u32;
    let mut m: u32;
    let head: u8;
    let mut hb: u8;

    if *half == 0 {
        let dix = *data.add(*di) >> 4;
        head = dix as u8;
    } else {
        let dix = *data.add(*di) & 0xf;
        head = dix as u8;
        *di += 1;
    }

    *half = 1-(*half);
    *res = 0;

    if head <= 8 {
        n = head as usize;
    } else {
        // leading ones, fill n half bytes in res
        n = (head - 8) as usize;
        mask = 0xf0000000;
        for i in 0..n {
            m = mask >> (4*i);
            *res = *res | m;
        }
    }

    if n == 8 {
        return Ok(());
    }

    if *di + ((8 - n) - (1 - *half)) / 2 >= max_di {
        return Err(From::from(ErrorKind::CorruptInputData));
    }

    for i in n..8 {
        if *half == 0 {
            let dix = *data.add(*di) >> 4;
            hb = dix as u8;
        } else {
            let dix = *data.add(*di) & 0xf;
            hb = dix as u8;
            *di += 1;
        }
        let hb32 = hb as u32;
        *res = *res | (hb32 << ((i-n)*4));
        *half = 1 - (*half);
    }

    Ok(())
}

// ENCODE/DECODE

/// Encodes double array using lossy conversion with value prediction.
///
/// The resulting binary is maximally 8 + size * 5 bytes, but much less if the
/// data is reasonably smooth on the first order.
///
/// This encoding is suitable for typical m/z or retention time binary arrays.
/// On a test set, the encoding was empirically show to be accurate to at
/// least 0.002 ppm.
///
/// Returns the number of encoded bytes.
///
/// * `data`        - Pointer to double array to be encoded
/// * `size`        - Size of double array
/// * `result`      - Pointer to resulting byte buffer
/// * `scaling`     - The scaling factor used for getting the fixed point repr.
pub unsafe extern "C" fn encode_linear(
    data: *const f64,
    data_size: usize,
    result: *mut u8,
    scaling: f64
)
    -> Result<usize>
{
    let mut ints: [i64; 3] = uninitialized();
    let mut extrapol: i64;

    encode_fixed_point(scaling, result);

    if data_size == 0 {
        return Ok(8);
    }

    ints[1] = (*data * scaling + 0.5) as i64;
    for i in 0..4 {
        let di = result.add(8+i);
        let xi = (ints[1] >> (i*8)) & 0xff;
        *di = xi as u8;
    }

    if data_size == 1 {
        return Ok(12);
    }

    ints[2] = (*data.add(1) * scaling + 0.5) as i64;
    for i in 0..4 {
        let di = result.add(12+i);
        let xi = (ints[2] >> (i*8)) & 0xff;
        *di = xi as u8;
    }

    let mut half_byte_count: usize = 0;
    let mut ri: usize = 16;
    let mut half_bytes: [u8; 10] = uninitialized();
    const I32_MIN: i64 = i32::min_value() as i64;
    const I32_MAX: i64 = i32::max_value() as i64;
    let mut diff: i32;

    for i in 2..data_size {
        ints[0] = ints[1];
        ints[1] = ints[2];
        if (*data.add(i) * scaling + 0.5) as i64 > i64::max_value() {
            return Err(From::from(ErrorKind::OverflowError));
        }

        ints[2] = ((*data.add(i)) * scaling + 0.5) as i64;
        extrapol = ints[1] + (ints[1] - ints[0]);

        if (ints[2] - extrapol > I32_MAX) || (ints[2] - extrapol < I32_MIN) {
            return Err(From::from(ErrorKind::OutOfRange));
        }

        diff = (ints[2] - extrapol) as i32;
        encode_int(diff as u32, &mut half_bytes[half_byte_count], &mut half_byte_count);

        for hbi in (1..half_byte_count).step_by(2) {
            let di = result.add(ri);
            let xi = (half_bytes[hbi-1] << 4) | (half_bytes[hbi] & 0xf);
            *di = xi as u8;
            ri += 1;
        }

        if half_byte_count % 2 != 0 {
            half_bytes[0] = half_bytes[half_byte_count-1];
            half_byte_count = 1;
        } else {
            half_byte_count = 0;
        }
    }

    if half_byte_count == 1 {
        let di = result.add(ri);
        let xi = half_bytes[0] << 4;
        *di = xi as u8;
        ri += 1;
    }

    Ok(ri)
}

/// Decodes data encoded by encode_linear.
///
/// The output size is guaranteed to be shorter or equal to
/// (|size| - 8) * 2.
///
/// Note that this method may throw an error if it deems the input data
/// to be corrupt, i.e. the last encoded int does not use the last byte
/// in the data. In addition the last encoded int need to use either the
/// last halfbyte, or the second last followed by a 0x0 halfbyte.
///
/// Returns the number of decoded doubles.
///
/// * `src`  - Pointer to bytes to be decoded.
/// * `size` - Size of byte array.
/// * `dst`  - Pointer to resulting double array.
pub unsafe extern "C" fn decode_linear(
    data: *const u8,
    data_size: usize,
    result: *mut f64
)
    -> Result<usize>
{
    // safety checks
    if data_size == 8 {
        return Ok(0);
    }

    if data_size < 8 {
        return Err(From::from(ErrorKind::CorruptInputData));
    }

    if data_size < 12 {
        return Err(From::from(ErrorKind::CorruptInputData));
    }

    let scaling = decode_fixed_point(data);
    let mut ints: [i64; 3] = uninitialized();
    let mut extrapol: i64;
    let mut init: i64;
    ints[1] = 0;
    for i in 0..4 {
        init = *data.add(8+i) as i64;
        let xi = (0xff & (init)) << (i*8);
        ints[1] = ints[1] | xi;
    }
    *result = (ints[1] as f64) / scaling;

    if data_size == 12 {
        return Ok(1);
    }

    if data_size < 16 {
        return Err(From::from(ErrorKind::CorruptInputData));
    }

    ints[2] = 0;
    for i in 0..4 {
        init = *data.add(12+i) as i64;
        let xi = (0xff & (init)) << (i*8);
        ints[2] = ints[2] | xi;
    }
    *result.add(1) = (ints[2] as f64) / scaling;

    let mut half: usize = 0;
    let mut ri: usize = 2;
    let mut di: usize = 16;
    let mut buff: u32 = uninitialized();
    let mut diff: i32;
    let mut y: i64;

    while di < data_size {
        if di == (data_size - 1) && half == 1 {
            if (*data.add(di) & 0xf) == 0x0 {
                break;
            }
        }

        ints[0] = ints[1];
        ints[1] = ints[2];
        decode_int(data, &mut di, data_size, &mut half, &mut buff)?;
        diff = buff as i32;

        extrapol = ints[1] + (ints[1] - ints[0]);
        y = extrapol + diff as i64;
        *result.add(ri) = (y as f64) / scaling;
        ri += 1;
        ints[2] = y;
    }

    Ok(ri)
}

// OPTIMAL

/// Macro for maximum value using partial ordering.
macro_rules! max {
    ($d0:ident, $d1:ident) => (if $d0 < $d1 { $d1 } else { $d0 })
}

/// Calculate the optimal scaling factor for Numpress compression.
pub unsafe extern "C" fn optimal_linear_scaling(
    data: *const f64,
    data_size: usize
)
    -> f64
{
    match data_size {
        0 => 0.,
        // 2147483647.0 == 0x7FFFFFFFl
        1 => (2147483647.0 / *data).floor(),
        _ => {
            let d0: f64 = *data;
            let d1: f64 = *data.add(1);
            let mut max_double: f64 = max!(d0, d1);
            let mut extrapol: f64;
            let mut diff: f64;

            for i in 2..data_size {
                let di: f64 = *data.add(i);
                let di_1: f64 = *data.add(i-1);
                let di_2: f64 = *data.add(i-2);
                extrapol = di_1 + (di_1 - di_2);
                diff = di - extrapol;
                let maxi = (diff.abs()+1.0).ceil();
                max_double = max!(max_double, maxi);
            }
            // 2147483647.0 == 0x7FFFFFFFl
            (2147483647.0 / max_double).floor()
        }
    }
}

}   // low_level

// API

/// Default scaling factor for compression.
pub const DEFAULT_SCALING: f64 = 10000.0;

/// High-level compressor for Numpress.
///
/// The recommended scaling factor is [`DEFAULT_SCALING`], and the optimal scaling
/// factor can be calculated via [`optimal_scaling`].
///
/// * `data`    - Slice of doubles to be encoded.
/// * `scaling` - Scaling factor used for getting the fixed point representation.
///
/// [`DEFAULT_SCALING`]: constant.DEFAULT_SCALING.html
/// [`optimal_scaling`]: fn.optimal_scaling.html
pub fn numpress_compress(data: &[f64], scaling: f64)
    -> Result<Vec<u8>>
{
    let mut vec: Vec<u8> = Vec::with_capacity(data.len() * 5 + 8);

    unsafe {
        let src: *const f64 = data.as_ptr();
        let dst: *mut u8 = vec.as_mut_ptr();
        let length = low_level::encode_linear(src, data.len(), dst, scaling)?;
        vec.set_len(length);
        vec.shrink_to_fit();
    }

    Ok(vec)
}

/// High-level decompressor for Numpress.
///
/// * `data`    - Slice of encoded doubles as bytes.
pub fn numpress_decompress(data: &[u8])
    -> Result<Vec<f64>>
{
    let mut vec: Vec<f64> = Vec::with_capacity((data.len() -8) * 2);

    unsafe {
        let src: *const u8 = data.as_ptr();
        let dst: *mut f64 = vec.as_mut_ptr();
        let length = low_level::decode_linear(src, data.len(), dst)?;
        vec.set_len(length);
        vec.shrink_to_fit();
    }

    Ok(vec)
}

/// Calculate the optimal, most-compressed scaling factor for compression.
///
/// * `data`    - Slice of doubles to be encoded.
pub fn optimal_scaling(data: &[f64])
    -> f64
{
    unsafe {
        low_level::optimal_linear_scaling(data.as_ptr(), data.len())
    }
}


#[cfg(test)]
mod tests {
    use rand::{thread_rng, Rng};
    use rand::distributions::Uniform;
    use std::alloc::{alloc, dealloc, Layout};
    use std::mem::{align_of, size_of, uninitialized};
    use super::*;

    macro_rules! assert_approx_list_eq {
        ($a:expr, $b:expr) => {
            assert_eq!($a.len(), $b.len());
            let mut iter = $a.iter().zip($b.iter());
            for (x, y) in iter {
                assert_approx_eq!(x, y);
            }
        };
        ($a:expr, $b:expr, $eps:expr) => {
            assert_eq!($a.len(), $b.len());
            let iter = $a.iter().zip($b.iter());
            for (x, y) in iter {
                assert_approx_eq!(x, y, $eps);
            }
        };
    }

    // FIXED POINT

    #[test]
    fn fixed_point_test() {
        unsafe {
            let x: f64 = 32.5;
            let mut xi: [u8; 8] = uninitialized();
            low_level::encode_fixed_point(x, xi.as_mut_ptr());
            assert_eq!(low_level::decode_fixed_point(xi.as_ptr()), x);

            let y: f64 = 1.2e-64;
            let mut yi: [u8; 8] = uninitialized();
            low_level::encode_fixed_point(y, yi.as_mut_ptr());
            assert_eq!(low_level::decode_fixed_point(yi.as_ptr()), y);
        }
    }

    #[test]
    fn fixed_point_heap_test() {
        unsafe {
            const SIZE: usize = size_of::<u8>() * 8;
            const ALIGN: usize = align_of::<u8>();
            let layout = Layout::from_size_align_unchecked(SIZE, ALIGN);

            let x: f64 = 32.5;
            let xi: *mut u8 = alloc(layout);
            low_level::encode_fixed_point(x, xi);
            assert_eq!(low_level::decode_fixed_point(xi), x);

            dealloc(xi, layout);
        }
    }

    // API

    #[test]
    fn compress_test() {
        // Check value compression with default scaling.
        let decoded: Vec<f64> = vec![100., 101., 102., 103.];
        let encoded: Vec<u8> = vec![64, 195, 136, 0, 0, 0, 0, 0, 64, 66, 15, 0, 80, 105, 15, 0, 136];
        let result = numpress_compress(&decoded, DEFAULT_SCALING).unwrap();
        assert_eq!(result, encoded);

        // Check value compression with optimal scaling.
        let encoded: Vec<u8> = vec![65, 116, 70, 248, 96, 0, 0, 0, 88, 144, 187, 126, 222, 255, 255, 127, 136];
        let result = numpress_compress(&decoded, 21262214.0).unwrap();
        assert_eq!(result, encoded);

        // Bug fix with custom input.
        let decoded: Vec<f64> = vec![472.36640759869624, 8161.255047730418, 31419.174861096908, 31340.37083086082, 11031.961448006856, 35019.3535619803, 22837.824611949254, 2076.226408785704, 23277.55357717535, 37604.579217858874, 34185.89109314591, 5077.6548386088325];
        let encoded: Vec<u8> = vec![64, 195, 136, 0, 0, 0, 0, 0, 208, 19, 72, 0, 6, 79, 221, 4, 25, 69, 167, 73, 152, 57, 23, 18, 155, 5, 49, 243, 0, 192, 7, 106, 16, 72, 240, 23, 174, 156, 9, 194, 234, 6, 200, 3, 9, 25, 137, 1, 126, 185, 240, 131, 198, 89, 96, 97, 11, 0];
        let result = numpress_compress(&decoded, DEFAULT_SCALING).unwrap();
        assert_eq!(result, encoded);
    }

    #[test]
    fn decompress_test() {
        // Check value decompression.
        let encoded: [u8; 17] = [64, 195, 136, 0, 0, 0, 0, 0, 64, 66, 15, 0, 80, 105, 15, 0, 136];
        let decoded: Vec<f64> = vec![100., 101., 102., 103.];
        let result = numpress_decompress(&encoded).unwrap();
        assert_eq!(result, decoded);

        // Bug fix with custom input.
        let encoded: Vec<u8> = vec![64, 195, 136, 0, 0, 0, 0, 0, 208, 19, 72, 0, 6, 79, 221, 4, 25, 69, 167, 73, 152, 57, 23, 18, 155, 5, 49, 243, 0, 192, 7, 106, 16, 72, 240, 23, 174, 156, 9, 194, 234, 6, 200, 3, 9, 25, 137, 1, 126, 185, 240, 131, 198, 89, 96, 97, 11, 0];
        let decoded: Vec<f64> = vec![472.3664, 8161.255, 31419.1749, 31340.3708, 11031.9614, 35019.3536, 22837.8246, 2076.2264, 23277.5536, 37604.5792, 34185.8911, 5077.6548];
        let result = numpress_decompress(&encoded).unwrap();
        assert_approx_list_eq!(result, decoded, 0.001);
    }

    #[test]
    fn optimal_scaling_test() {
        // Check optimal fixed point detection
        let decoded: [f64; 4] = [100., 101., 102., 103.];
        assert_eq!(optimal_scaling(&decoded), 21262214.0);
    }

    #[test]
    #[ignore]
    fn fuzz_test() {
        // fuzz with random integers to ensure minimal loss and no memory corruption
        let mut rng = thread_rng();
        let dist = Uniform::new(0f64, 55000f64);
        for _ in 0..10000 {
            let length: usize = rng.gen_range(0, 10000);
            let input: Vec<f64> = rng.sample_iter(&dist).take(length).collect();
            let encoded = numpress_compress(input.as_slice(), DEFAULT_SCALING).unwrap();
            let decoded = numpress_decompress(encoded.as_slice()).unwrap();
            assert_approx_list_eq!(input, decoded, 0.0001);
        }
    }
}
