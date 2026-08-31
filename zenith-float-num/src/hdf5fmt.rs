//! HDF5 array I/O via crates.io `hdf5-rust`. No libhdf5.
//!
//! IEEE arrays store native binary32 / binary64 bit lanes. `ExactNumArray`
//! stores one opaque cell per entry: a big-endian `u32` payload length, then
//! [`ExactNum::to_bytes`], zero-padded to a shared element size.
//!
//! CSV is unchanged. Requires features `std` and `hdf5`.

use crate::defs::DEFAULT_P;
use crate::ieee_soft::{ExactNumArray, Ieee32, Ieee32Array, Ieee64Array};
use crate::Error;
use crate::ExactNum;
use crate::Hdf5File;
use alloc::vec::Vec;

fn map_hdf5(_err: hdf5_rust::HDF5Error) -> Error {
    Error::InvalidArgument
}

fn dataspace(rows: usize, cols: usize) -> Result<Vec<usize>, Error> {
    if rows == 0 || cols == 0 {
        return Err(Error::InvalidArgument);
    }
    if rows == 1 {
        Ok(alloc::vec![cols])
    } else {
        Ok(alloc::vec![rows, cols])
    }
}

fn from_dataspace(shape: &[usize]) -> Result<(usize, usize), Error> {
    match *shape {
        [n] => Ok((1, n)),
        [r, c] => Ok((r, c)),
        _ => Err(Error::InvalidArgument),
    }
}

fn open_or_create(path: &std::path::Path) -> Result<Hdf5File, Error> {
    if path.exists() {
        Hdf5File::open(path).map_err(map_hdf5)
    } else {
        Ok(Hdf5File::create())
    }
}

impl Ieee64Array {
    /// Write this array into `file` at `dataset` (rank-1 when the shape is a
    /// row vector `(1, n)`, otherwise rank-2).
    pub fn write_hdf5(&self, file: &mut Hdf5File, dataset: &str) -> Result<(), Error> {
        let (rows, cols) = self.shape();
        let shape = dataspace(rows, cols)?;
        file.write_ieee64(dataset, &shape, self.as_bits())
            .map_err(map_hdf5)
    }

    /// Read an IEEE binary64 dataset from an open file.
    pub fn read_hdf5(file: &Hdf5File, dataset: &str) -> Result<Self, Error> {
        let (shape, bits) = file.read_ieee64(dataset).map_err(map_hdf5)?;
        let (rows, cols) = from_dataspace(&shape)?;
        Self::from_parts(rows, cols, bits)
    }

    /// Append this array's rows onto an existing IEEE binary64 dataset.
    pub fn append_hdf5_file(&self, file: &mut Hdf5File, dataset: &str) -> Result<(), Error> {
        let (rows, cols) = self.shape();
        let shape = dataspace(rows, cols)?;
        file.append_ieee64(dataset, &shape, self.as_bits())
            .map_err(map_hdf5)
    }

    /// Create or open `path` and write `dataset`.
    pub fn to_hdf5<P: AsRef<std::path::Path>>(&self, path: P, dataset: &str) -> Result<(), Error> {
        let path = path.as_ref();
        let mut file = open_or_create(path)?;
        self.write_hdf5(&mut file, dataset)?;
        file.save(path).map_err(map_hdf5)
    }

    /// Read `dataset` from an HDF5 file at `path`.
    pub fn from_hdf5<P: AsRef<std::path::Path>>(path: P, dataset: &str) -> Result<Self, Error> {
        let file = Hdf5File::open(path).map_err(map_hdf5)?;
        Self::read_hdf5(&file, dataset)
    }

    /// Open `path` and append this array's rows onto `dataset`.
    pub fn append_hdf5<P: AsRef<std::path::Path>>(
        &self,
        path: P,
        dataset: &str,
    ) -> Result<(), Error> {
        let path = path.as_ref();
        let mut file = Hdf5File::open(path).map_err(map_hdf5)?;
        self.append_hdf5_file(&mut file, dataset)?;
        file.save(path).map_err(map_hdf5)
    }
}

impl Ieee32Array {
    /// Write this array into `file` at `dataset`.
    pub fn write_hdf5(&self, file: &mut Hdf5File, dataset: &str) -> Result<(), Error> {
        let (rows, cols) = self.shape();
        let shape = dataspace(rows, cols)?;
        file.write_ieee32(dataset, &shape, self.as_bits())
            .map_err(map_hdf5)
    }

    /// Read an IEEE binary32 dataset from an open file.
    pub fn read_hdf5(file: &Hdf5File, dataset: &str) -> Result<Self, Error> {
        let (shape, bits) = file.read_ieee32(dataset).map_err(map_hdf5)?;
        let (rows, cols) = from_dataspace(&shape)?;
        let mut vals = Vec::new();
        vals.try_reserve_exact(bits.len())?;
        for b in bits {
            vals.push(Ieee32::from_bits(b));
        }
        Self::from_shape(rows, cols, &vals).ok_or(Error::InvalidArgument)
    }

    /// Create or open `path` and write `dataset`.
    pub fn to_hdf5<P: AsRef<std::path::Path>>(&self, path: P, dataset: &str) -> Result<(), Error> {
        let path = path.as_ref();
        let mut file = open_or_create(path)?;
        self.write_hdf5(&mut file, dataset)?;
        file.save(path).map_err(map_hdf5)
    }

    /// Read `dataset` from an HDF5 file at `path`.
    pub fn from_hdf5<P: AsRef<std::path::Path>>(path: P, dataset: &str) -> Result<Self, Error> {
        let file = Hdf5File::open(path).map_err(map_hdf5)?;
        Self::read_hdf5(&file, dataset)
    }
}

impl ExactNumArray {
    /// Pack each cell as `u32` BE length + [`ExactNum::to_bytes`], padded.
    pub fn write_hdf5(&self, file: &mut Hdf5File, dataset: &str) -> Result<(), Error> {
        let (rows, cols) = self.shape();
        let shape = dataspace(rows, cols)?;
        let n = self.len();
        let mut payloads = Vec::new();
        payloads.try_reserve_exact(n)?;
        let mut max_len = 0usize;
        for i in 0..n {
            let v = self.get(i).ok_or(Error::InvalidArgument)?;
            let bytes = v.to_bytes()?;
            if bytes.len() > max_len {
                max_len = bytes.len();
            }
            payloads.push(bytes);
        }
        let prefix = 4usize;
        let elem = prefix.checked_add(max_len).ok_or(Error::InvalidArgument)?;
        let nbytes = n.checked_mul(elem).ok_or(Error::InvalidArgument)?;
        let mut data = Vec::new();
        data.try_reserve_exact(nbytes)?;
        data.resize(nbytes, 0);
        for (i, payload) in payloads.iter().enumerate() {
            let plen = u32::try_from(payload.len()).map_err(|_| Error::InvalidArgument)?;
            let off = i.checked_mul(elem).ok_or(Error::InvalidArgument)?;
            data[off..off + 4].copy_from_slice(&plen.to_be_bytes());
            data[off + 4..off + 4 + payload.len()].copy_from_slice(payload);
        }
        file.write_opaque(dataset, &shape, elem, &data)
            .map_err(map_hdf5)
    }

    /// Decode an opaque dataset of length-prefixed [`ExactNum`] records.
    pub fn read_hdf5(file: &Hdf5File, dataset: &str) -> Result<Self, Error> {
        let (shape, elem, data) = file.read_opaque(dataset).map_err(map_hdf5)?;
        let (rows, cols) = from_dataspace(&shape)?;
        if elem < 4 {
            return Err(Error::InvalidArgument);
        }
        let n = rows.checked_mul(cols).ok_or(Error::InvalidArgument)?;
        let need = n.checked_mul(elem).ok_or(Error::InvalidArgument)?;
        if data.len() != need {
            return Err(Error::InvalidArgument);
        }
        let mut vals = Vec::new();
        vals.try_reserve_exact(n)?;
        for i in 0..n {
            let off = i.checked_mul(elem).ok_or(Error::InvalidArgument)?;
            let cell = &data[off..off + elem];
            let plen = u32::from_be_bytes([cell[0], cell[1], cell[2], cell[3]]) as usize;
            let end = 4usize.checked_add(plen).ok_or(Error::InvalidArgument)?;
            if end > elem {
                return Err(Error::InvalidArgument);
            }
            vals.push(ExactNum::from_bytes(&cell[4..end])?);
        }
        let p = vals.iter().find_map(|v| v.precision()).unwrap_or(DEFAULT_P);
        Self::from_parts(p, rows, cols, vals)
    }

    /// Create or open `path` and write `dataset`.
    pub fn to_hdf5<P: AsRef<std::path::Path>>(&self, path: P, dataset: &str) -> Result<(), Error> {
        let path = path.as_ref();
        let mut file = open_or_create(path)?;
        self.write_hdf5(&mut file, dataset)?;
        file.save(path).map_err(map_hdf5)
    }

    /// Read `dataset` from an HDF5 file at `path`.
    pub fn from_hdf5<P: AsRef<std::path::Path>>(path: P, dataset: &str) -> Result<Self, Error> {
        let file = Hdf5File::open(path).map_err(map_hdf5)?;
        Self::read_hdf5(&file, dataset)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ieee_soft::Ieee64;
    use crate::ExactNum;
    use crate::Hdf5File;
    use crate::Word;

    fn tmp(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("zenith_hdf5_{}_{}.h5", std::process::id(), name))
    }

    fn cleanup(path: &std::path::Path) {
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn gold_ieee64_100x3_bit_identical() {
        let mut vals = Vec::with_capacity(300);
        for i in 0..300i32 {
            vals.push(Ieee64::from_i32(i));
        }
        let a = Ieee64Array::from_shape(100, 3, &vals).unwrap();
        let path = tmp("ieee64_100x3");
        cleanup(&path);
        a.to_hdf5(&path, "data").unwrap();
        let b = Ieee64Array::from_hdf5(&path, "data").unwrap();
        cleanup(&path);
        assert_eq!(b.shape(), (100, 3));
        assert_eq!(a.as_bits(), b.as_bits());
        assert_eq!(b.as_bits().len(), 300);
    }

    #[test]
    fn gold_exact_50x50_p256() {
        const P: usize = 256;
        let n = 50 * 50;
        let mut vals = Vec::with_capacity(n);
        for i in 0..n {
            vals.push(ExactNum::from_word((i as Word).wrapping_add(1), P));
        }
        let a = ExactNumArray::from_shape(P, 50, 50, &vals).unwrap();
        let path = tmp("exact_50x50");
        cleanup(&path);
        a.to_hdf5(&path, "exact").unwrap();
        let b = ExactNumArray::from_hdf5(&path, "exact").unwrap();
        cleanup(&path);
        assert_eq!(b.shape(), (50, 50));
        assert_eq!(b.precision(), P);
        assert_eq!(b.len(), n);
        for i in 0..n {
            assert_eq!(
                a.get(i).unwrap().cmp(b.get(i).unwrap()),
                Some(0),
                "cell {i}"
            );
        }
    }

    #[test]
    fn gold_ieee32_1d_1000_bit_identical() {
        let bits: Vec<u32> = (0..1000u32)
            .map(|i| 0x3F80_0000u32.wrapping_add(i))
            .collect();
        let a = Ieee32Array::from_bits(&bits);
        assert_eq!(a.shape(), (1, 1000));
        let path = tmp("ieee32_1d_1000");
        cleanup(&path);
        a.to_hdf5(&path, "vec").unwrap();
        let b = Ieee32Array::from_hdf5(&path, "vec").unwrap();
        cleanup(&path);
        assert_eq!(b.shape(), (1, 1000));
        assert_eq!(a.as_bits(), b.as_bits());
    }

    #[test]
    fn gold_nested_group_results_data() {
        let a = Ieee64Array::from_shape(
            2,
            3,
            &[
                Ieee64::from_i32(1),
                Ieee64::from_i32(2),
                Ieee64::from_i32(3),
                Ieee64::from_i32(4),
                Ieee64::from_i32(5),
                Ieee64::from_i32(6),
            ],
        )
        .unwrap();
        let path = tmp("nested_results");
        cleanup(&path);
        let mut file = Hdf5File::create();
        file.create_group("results").unwrap();
        a.write_hdf5(&mut file, "results/data").unwrap();
        file.save(&path).unwrap();
        let b = Ieee64Array::from_hdf5(&path, "results/data").unwrap();
        cleanup(&path);
        assert_eq!(a.as_bits(), b.as_bits());
        assert_eq!(b.shape(), (2, 3));
    }

    #[test]
    fn gold_append_10x3_to_20x3() {
        let mut first = Vec::with_capacity(30);
        let mut second = Vec::with_capacity(30);
        for i in 0..30i32 {
            first.push(Ieee64::from_i32(i));
            second.push(Ieee64::from_i32(i + 30));
        }
        let a = Ieee64Array::from_shape(10, 3, &first).unwrap();
        let b = Ieee64Array::from_shape(10, 3, &second).unwrap();
        let path = tmp("append_10x3");
        cleanup(&path);
        a.to_hdf5(&path, "data").unwrap();
        b.append_hdf5(&path, "data").unwrap();
        let c = Ieee64Array::from_hdf5(&path, "data").unwrap();
        cleanup(&path);
        assert_eq!(c.shape(), (20, 3));
        assert_eq!(&c.as_bits()[..30], a.as_bits());
        assert_eq!(&c.as_bits()[30..], b.as_bits());
    }

    #[test]
    fn gold_wrong_dataset_name_is_err() {
        let a = Ieee64Array::from_values(&[Ieee64::from_i32(1)]);
        let path = tmp("wrong_name");
        cleanup(&path);
        a.to_hdf5(&path, "data").unwrap();
        assert!(Ieee64Array::from_hdf5(&path, "missing").is_err());
        assert!(ExactNumArray::from_hdf5(&path, "missing").is_err());
        cleanup(&path);
    }
}
