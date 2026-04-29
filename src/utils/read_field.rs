//! Simple wrappers for BAM auxiliary tags.
//!
//! This module provides owned representations of selected `rust-htslib`
//! auxiliary tag values so they can be stored independently of the source BAM
//! record. Array-valued auxiliary tags are intentionally rejected.

use crate::Error;
use rust_htslib::bam::record::Aux;
use serde::de::{self, MapAccess, Visitor};
use serde::ser::SerializeMap as _;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::borrow::Borrow;
use std::collections::BTreeMap;
use std::collections::btree_map::{Entry, Iter, IterMut, Keys, Values, ValuesMut};
use std::fmt;
use std::ops::Index;

/// Owned representation of simple BAM auxiliary tag values.
///
/// `SimpleAux` is intentionally limited to scalar BAM aux values that map
/// cleanly to common serialized forms.
///
/// # Examples
///
/// ```
/// use nanalogue_core::SimpleAux;
///
/// assert_eq!(serde_json::to_string(&SimpleAux::I64(1))?, "1");
/// assert_eq!(serde_json::to_string(&SimpleAux::Double(1.25))?, "1.25");
/// assert_eq!(
///     serde_json::to_string(&SimpleAux::String(String::from("x")))?,
///     r#""x""#
/// );
/// # Ok::<(), serde_json::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[expect(
    clippy::exhaustive_enums,
    reason = "public enum is intended for direct matching in downstream code"
)]
#[serde(untagged)]
pub enum SimpleAux {
    /// Signed 64-bit integer.
    I64(i64),
    /// 64-bit floating point value.
    Double(f64),
    /// UTF-8 string value.
    String(String),
}

impl From<i64> for SimpleAux {
    fn from(value: i64) -> Self {
        Self::I64(value)
    }
}

impl From<f64> for SimpleAux {
    fn from(value: f64) -> Self {
        Self::Double(value)
    }
}

impl From<String> for SimpleAux {
    fn from(value: String) -> Self {
        Self::String(value)
    }
}

impl From<&str> for SimpleAux {
    fn from(value: &str) -> Self {
        Self::String(value.to_owned())
    }
}

/// Formats the contained scalar value without adding JSON quoting.
///
/// # Examples
///
/// ```
/// use nanalogue_core::SimpleAux;
///
/// assert_eq!(SimpleAux::I64(1).to_string(), "1");
/// assert_eq!(SimpleAux::Double(1.25).to_string(), "1.25");
/// assert_eq!(SimpleAux::String(String::from("x")).to_string(), "x");
/// ```
impl fmt::Display for SimpleAux {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.clone() {
            Self::I64(value) => write!(f, "{value}"),
            Self::Double(value) => write!(f, "{value}"),
            Self::String(value) => write!(f, "{value}"),
        }
    }
}

/// Converts supported `rust-htslib` BAM aux values into owned `SimpleAux`
/// values.
///
/// Conversion policy:
/// - `Aux::Char` becomes `SimpleAux::String`
/// - integer-valued aux variants become `SimpleAux::I64`
/// - floating-point aux variants become `SimpleAux::Double`
/// - `Aux::String` becomes `SimpleAux::String`
/// - array-valued aux tags are rejected
/// - `Aux::HexByteArray` is rejected
///
/// # Examples
///
/// ```
/// use nanalogue_core::{Error, SimpleAux};
/// use rust_htslib::bam::record::Aux;
///
/// assert_eq!(SimpleAux::try_from(Aux::Char(b'A'))?, SimpleAux::String(String::from("A")));
/// assert_eq!(SimpleAux::try_from(Aux::I32(-7))?, SimpleAux::I64(-7));
/// assert_eq!(SimpleAux::try_from(Aux::Double(1.5))?, SimpleAux::Double(1.5));
/// assert!(matches!(
///     SimpleAux::try_from(Aux::HexByteArray("CAFE")),
///     Err(Error::NotImplemented(_))
/// ));
/// # Ok::<(), Error>(())
/// ```
impl TryFrom<Aux<'_>> for SimpleAux {
    type Error = Error;

    fn try_from(value: Aux<'_>) -> Result<Self, Self::Error> {
        match value {
            Aux::Char(v) => Ok(Self::String(char::from(v).to_string())),
            Aux::I8(v) => Ok(Self::I64(i64::from(v))),
            Aux::U8(v) => Ok(Self::I64(i64::from(v))),
            Aux::I16(v) => Ok(Self::I64(i64::from(v))),
            Aux::U16(v) => Ok(Self::I64(i64::from(v))),
            Aux::I32(v) => Ok(Self::I64(i64::from(v))),
            Aux::U32(v) => Ok(Self::I64(i64::from(v))),
            Aux::Float(v) => Ok(Self::Double(f64::from(v))),
            Aux::Double(v) => Ok(Self::Double(v)),
            Aux::String(v) => Ok(Self::String(v.to_owned())),
            Aux::HexByteArray(_) => Err(Error::NotImplemented(
                "hex byte array BAM aux tags are not supported in `SimpleAux`".to_owned(),
            )),
            Aux::ArrayI8(_)
            | Aux::ArrayU8(_)
            | Aux::ArrayI16(_)
            | Aux::ArrayU16(_)
            | Aux::ArrayI32(_)
            | Aux::ArrayU32(_)
            | Aux::ArrayFloat(_) => Err(Error::NotImplemented(
                "array-valued BAM aux tags are not supported in `SimpleAux`".to_owned(),
            )),
        }
    }
}

/// Owned map of BAM auxiliary tags keyed by their two-byte names.
///
/// `ReadField` wraps a `BTreeMap<[u8; 2], SimpleAux>` and preserves sorted key
/// order.
///
/// # Examples
///
/// ```
/// use nanalogue_core::{ReadField, SimpleAux};
///
/// let mut read_field = ReadField::new();
/// assert_eq!(read_field.insert(*b"ts", 1i64), None);
/// assert_eq!(read_field.insert(*b"ul", "hello"), None);
///
/// assert_eq!(read_field.get(&b"ts"[..]), Some(&SimpleAux::I64(1)));
/// assert_eq!(read_field[&b"ul"[..]], SimpleAux::String(String::from("hello")));
/// assert_eq!(serde_json::to_string(&read_field)?, r#"{"ts":1,"ul":"hello"}"#);
/// # Ok::<(), serde_json::Error>(())
/// ```
///
/// Additional ergonomic construction and bulk-update patterns:
///
/// ```
/// use nanalogue_core::{ReadField, SimpleAux};
///
/// let from_array = ReadField::from([(*b"ts", 1i64), (*b"ul", 2i64)]);
/// let from_iter = [(*b"ts", 1i64), (*b"ul", 2i64)]
///     .into_iter()
///     .collect::<ReadField>();
///
/// let mut extended = ReadField::new();
/// extended.extend([(*b"ts", 1i64), (*b"ul", 2i64)]);
///
/// let borrowed_entries = [
///     (*b"ts", SimpleAux::I64(1)),
///     (*b"ul", SimpleAux::I64(2)),
/// ];
/// let mut extended_by_ref = ReadField::new();
/// extended_by_ref.extend(borrowed_entries.iter().map(|entry| (&entry.0, &entry.1)));
///
/// assert_eq!(from_array, from_iter);
/// assert_eq!(from_iter, extended);
/// assert_eq!(extended, extended_by_ref);
/// ```
///
/// Indexing is supported for borrowed keys, but like `BTreeMap` indexing it
/// will panic if the key is absent.
///
/// ```
/// use nanalogue_core::{ReadField, SimpleAux};
///
/// let read_field = ReadField::from([(*b"ts", 1i64)]);
/// assert_eq!(read_field[&b"ts"[..]], SimpleAux::I64(1));
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
#[expect(
    clippy::exhaustive_structs,
    reason = "newtype wrapper is intended for direct construction in downstream code"
)]
pub struct ReadField(pub BTreeMap<[u8; 2], SimpleAux>);

impl ReadField {
    /// Makes a new, empty `ReadField`.
    ///
    /// Does not allocate anything on its own.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::ReadField;
    ///
    /// let mut read_field = ReadField::new();
    ///
    /// assert_eq!(read_field.insert(*b"ts", 1i64), None);
    /// assert_eq!(read_field.len(), 1);
    /// ```
    #[must_use]
    pub const fn new() -> Self {
        Self(BTreeMap::new())
    }

    /// Insert a BAM auxiliary tag value, returning any previous value.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::{ReadField, SimpleAux};
    ///
    /// let mut read_field = ReadField::new();
    /// assert_eq!(read_field.insert(*b"ts", 1i64), None);
    /// assert_eq!(read_field.insert(*b"ul", "hello"), None);
    /// assert_eq!(read_field.insert(*b"ts", "updated"), Some(SimpleAux::I64(1)));
    ///
    /// assert_eq!(
    ///     read_field.get(&b"ul"[..]),
    ///     Some(&SimpleAux::String(String::from("hello")))
    /// );
    /// ```
    #[must_use]
    pub fn insert<U>(&mut self, key: [u8; 2], value: U) -> Option<SimpleAux>
    where
        U: Into<SimpleAux>,
    {
        self.0.insert(key, value.into())
    }

    /// Clears the map, removing all elements.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::ReadField;
    ///
    /// let mut read_field = ReadField::default();
    /// let _previous = read_field.insert(*b"ts", 1i64);
    /// read_field.clear();
    /// assert!(read_field.is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Returns the number of elements in the map.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::ReadField;
    ///
    /// let mut read_field = ReadField::default();
    /// assert_eq!(read_field.len(), 0);
    /// let _previous = read_field.insert(*b"ts", 1i64);
    /// assert_eq!(read_field.len(), 1);
    /// ```
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns `true` if the map contains no elements.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::ReadField;
    ///
    /// let mut read_field = ReadField::default();
    /// assert!(read_field.is_empty());
    /// let _previous = read_field.insert(*b"ts", 1i64);
    /// assert!(!read_field.is_empty());
    /// ```
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns `true` if the map contains a value for the specified key.
    ///
    /// The key may be any borrowed form of the map's key type, but the ordering
    /// on the borrowed form must match the ordering on the key type.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::ReadField;
    ///
    /// let mut read_field = ReadField::default();
    /// let _previous = read_field.insert(*b"ts", 1i64);
    /// assert!(read_field.contains_key(&b"ts"[..]));
    /// assert!(!read_field.contains_key(&b"ul"[..]));
    /// ```
    #[must_use]
    pub fn contains_key<Q>(&self, key: &Q) -> bool
    where
        [u8; 2]: Borrow<Q> + Ord,
        Q: Ord + ?Sized,
    {
        self.0.contains_key(key)
    }

    /// Returns a reference to the value corresponding to the key.
    ///
    /// The key may be any borrowed form of the map's key type, but the ordering
    /// on the borrowed form must match the ordering on the key type.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::{ReadField, SimpleAux};
    ///
    /// let mut read_field = ReadField::new();
    /// assert_eq!(read_field.insert(*b"ts", 1i64), None);
    /// assert_eq!(read_field.get(&b"ts"[..]), Some(&SimpleAux::I64(1)));
    /// assert_eq!(read_field.get(&b"ul"[..]), None);
    /// ```
    #[must_use]
    pub fn get<Q>(&self, key: &Q) -> Option<&SimpleAux>
    where
        [u8; 2]: Borrow<Q> + Ord,
        Q: Ord + ?Sized,
    {
        self.0.get(key)
    }

    /// Returns the key-value pair corresponding to the supplied key.
    ///
    /// This can be useful for getting a reference to the stored key value from
    /// a borrowed lookup key, or for getting a key reference with the same
    /// lifetime as the collection.
    ///
    /// The supplied key may be any borrowed form of the map's key type, but the
    /// ordering on the borrowed form must match the ordering on the key type.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::{ReadField, SimpleAux};
    ///
    /// let mut read_field = ReadField::new();
    /// assert_eq!(read_field.insert(*b"ts", 1i64), None);
    /// assert_eq!(
    ///     read_field.get_key_value(&b"ts"[..]),
    ///     Some((&*b"ts", &SimpleAux::I64(1))),
    /// );
    /// assert_eq!(read_field.get_key_value(&b"ul"[..]), None);
    /// ```
    #[must_use]
    pub fn get_key_value<Q>(&self, key: &Q) -> Option<(&[u8; 2], &SimpleAux)>
    where
        [u8; 2]: Borrow<Q> + Ord,
        Q: Ord + ?Sized,
    {
        self.0.get_key_value(key)
    }

    /// Returns a mutable reference to the value corresponding to the key.
    ///
    /// The key may be any borrowed form of the map's key type, but the ordering
    /// on the borrowed form must match the ordering on the key type.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::{ReadField, SimpleAux};
    ///
    /// let mut read_field = ReadField::new();
    /// assert_eq!(read_field.insert(*b"ts", 1i64), None);
    /// if let Some(value) = read_field.get_mut(&b"ts"[..]) {
    ///     *value = SimpleAux::I64(2);
    /// }
    /// assert_eq!(read_field[&b"ts"[..]], SimpleAux::I64(2));
    /// ```
    pub fn get_mut<Q>(&mut self, key: &Q) -> Option<&mut SimpleAux>
    where
        [u8; 2]: Borrow<Q> + Ord,
        Q: Ord + ?Sized,
    {
        self.0.get_mut(key)
    }

    /// Gets an iterator over the entries of the map, sorted by key.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::{ReadField, SimpleAux};
    ///
    /// let mut read_field = ReadField::new();
    /// assert_eq!(read_field.insert(*b"cc", "c"), None);
    /// assert_eq!(read_field.insert(*b"bb", "b"), None);
    /// assert_eq!(read_field.insert(*b"aa", "a"), None);
    ///
    /// let first = read_field.iter().next();
    /// assert_eq!(first, Some((&*b"aa", &SimpleAux::String(String::from("a")))));
    /// ```
    pub fn iter(&self) -> Iter<'_, [u8; 2], SimpleAux> {
        self.0.iter()
    }

    /// Gets a mutable iterator over the entries of the map, sorted by key.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::{ReadField, SimpleAux};
    ///
    /// let mut read_field = ReadField::from([(*b"aa", 1i64), (*b"bb", 2i64), (*b"cc", 3i64)]);
    ///
    /// for (key, value) in read_field.iter_mut() {
    ///     if key != b"aa" {
    ///         *value = match value.clone() {
    ///             SimpleAux::I64(n) => SimpleAux::I64(n + 10),
    ///             other => other,
    ///         };
    ///     }
    /// }
    ///
    /// assert_eq!(read_field[&b"aa"[..]], SimpleAux::I64(1));
    /// assert_eq!(read_field[&b"bb"[..]], SimpleAux::I64(12));
    /// assert_eq!(read_field[&b"cc"[..]], SimpleAux::I64(13));
    /// ```
    pub fn iter_mut(&mut self) -> IterMut<'_, [u8; 2], SimpleAux> {
        self.0.iter_mut()
    }

    /// Gets an iterator over the keys of the map, in sorted order.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::ReadField;
    ///
    /// let mut read_field = ReadField::new();
    /// assert_eq!(read_field.insert(*b"bb", "b"), None);
    /// assert_eq!(read_field.insert(*b"aa", "a"), None);
    ///
    /// let keys: Vec<[u8; 2]> = read_field.keys().copied().collect();
    /// assert_eq!(keys, [*b"aa", *b"bb"]);
    /// ```
    pub fn keys(&self) -> Keys<'_, [u8; 2], SimpleAux> {
        self.0.keys()
    }

    /// Gets the given key's corresponding entry in the map for in-place
    /// manipulation.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::{ReadField, SimpleAux};
    ///
    /// let mut counts = ReadField::default();
    ///
    /// for tag in [*b"ts", *b"ul", *b"ts"] {
    ///     counts
    ///         .entry(tag)
    ///         .and_modify(|value| {
    ///             if let SimpleAux::I64(curr) = value {
    ///                 *curr += 1;
    ///             }
    ///         })
    ///         .or_insert(SimpleAux::I64(1));
    /// }
    ///
    /// assert_eq!(counts[&b"ts"[..]], SimpleAux::I64(2));
    /// assert_eq!(counts[&b"ul"[..]], SimpleAux::I64(1));
    /// ```
    pub fn entry(&mut self, key: [u8; 2]) -> Entry<'_, [u8; 2], SimpleAux> {
        self.0.entry(key)
    }

    /// Gets an iterator over the values of the map, in order by key.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::{ReadField, SimpleAux};
    ///
    /// let mut read_field = ReadField::default();
    /// let _previous = read_field.insert(*b"ts", "hello");
    /// let _previous = read_field.insert(*b"ul", "goodbye");
    ///
    /// let values: Vec<SimpleAux> = read_field.values().cloned().collect();
    /// assert_eq!(values, [SimpleAux::String(String::from("hello")), SimpleAux::String(String::from("goodbye"))]);
    /// ```
    pub fn values(&self) -> Values<'_, [u8; 2], SimpleAux> {
        self.0.values()
    }

    /// Gets a mutable iterator over the values of the map, in order by key.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::{ReadField, SimpleAux};
    ///
    /// let mut read_field = ReadField::default();
    /// let _previous = read_field.insert(*b"ts", String::from("hello"));
    /// let _previous = read_field.insert(*b"ul", String::from("goodbye"));
    ///
    /// for value in read_field.values_mut() {
    ///     if let SimpleAux::String(text) = value {
    ///         text.push('!');
    ///     }
    /// }
    ///
    /// let values: Vec<SimpleAux> = read_field.values().cloned().collect();
    /// assert_eq!(
    ///     values,
    ///     [
    ///         SimpleAux::String(String::from("hello!")),
    ///         SimpleAux::String(String::from("goodbye!")),
    ///     ]
    /// );
    /// ```
    pub fn values_mut(&mut self) -> ValuesMut<'_, [u8; 2], SimpleAux> {
        self.0.values_mut()
    }
}

/// Extends a `ReadField` from borrowed key-value pairs by copying the two-byte
/// tag and cloning the `SimpleAux` value.
impl<'a> Extend<(&'a [u8; 2], &'a SimpleAux)> for ReadField {
    fn extend<T: IntoIterator<Item = (&'a [u8; 2], &'a SimpleAux)>>(&mut self, iter: T) {
        self.0
            .extend(iter.into_iter().map(|(&key, value)| (key, value.clone())));
    }
}

/// Extends a `ReadField` from owned key-value pairs, converting each value via
/// `Into<SimpleAux>`.
impl<U> Extend<([u8; 2], U)> for ReadField
where
    U: Into<SimpleAux>,
{
    fn extend<T: IntoIterator<Item = ([u8; 2], U)>>(&mut self, iter: T) {
        self.0
            .extend(iter.into_iter().map(|(key, value)| (key, value.into())));
    }
}

/// Builds a `ReadField` from an array of owned key-value pairs, converting
/// values via `Into<SimpleAux>`.
impl<U, const N: usize> From<[([u8; 2], U); N]> for ReadField
where
    U: Into<SimpleAux>,
{
    fn from(entries: [([u8; 2], U); N]) -> Self {
        Self(
            entries
                .into_iter()
                .map(|(key, value)| (key, value.into()))
                .collect(),
        )
    }
}

/// Collects owned key-value pairs into a `ReadField`, converting values via
/// `Into<SimpleAux>`.
impl<U> FromIterator<([u8; 2], U)> for ReadField
where
    U: Into<SimpleAux>,
{
    fn from_iter<T: IntoIterator<Item = ([u8; 2], U)>>(iter: T) -> Self {
        Self(
            iter.into_iter()
                .map(|(key, value)| (key, value.into()))
                .collect(),
        )
    }
}

impl<'a> IntoIterator for &'a ReadField {
    type Item = (&'a [u8; 2], &'a SimpleAux);
    type IntoIter = Iter<'a, [u8; 2], SimpleAux>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a> IntoIterator for &'a mut ReadField {
    type Item = (&'a [u8; 2], &'a mut SimpleAux);
    type IntoIter = IterMut<'a, [u8; 2], SimpleAux>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

/// Supports `read_field[&key]` lookups for borrowed key forms.
///
/// Like `BTreeMap` indexing, this panics if the key is not present.
impl<Q> Index<&Q> for ReadField
where
    [u8; 2]: Borrow<Q> + Ord,
    Q: Ord + ?Sized,
{
    type Output = SimpleAux;

    fn index(&self, index: &Q) -> &Self::Output {
        self.0
            .get(index)
            .expect("ReadField index key is expected to exist")
    }
}

impl Serialize for ReadField {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut map = serializer.serialize_map(Some(self.0.len()))?;
        for (tag, value) in &self.0 {
            let key = std::str::from_utf8(tag).map_err(serde::ser::Error::custom)?;
            map.serialize_entry(key, value)?;
        }
        map.end()
    }
}

impl<'de> Deserialize<'de> for ReadField {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct ReadFieldVisitor;

        impl<'de> Visitor<'de> for ReadFieldVisitor {
            type Value = ReadField;

            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str("a JSON object with two-character BAM aux tag keys")
            }

            fn visit_map<A>(self, mut map: A) -> Result<Self::Value, A::Error>
            where
                A: MapAccess<'de>,
            {
                let mut output = ReadField::default();
                while let Some((key, value)) = map.next_entry::<String, SimpleAux>()? {
                    let key_bytes = key.as_bytes();
                    let tag_name: [u8; 2] = key_bytes.try_into().map_err(|_error| {
                        de::Error::custom(format!(
                            "BAM aux tag keys must be exactly 2 bytes, got `{key}`"
                        ))
                    })?;
                    if output.insert(tag_name, value).is_some() {
                        return Err(de::Error::custom(format!(
                            "duplicate BAM aux tag key encountered during deserialization: `{key}`"
                        )));
                    }
                }
                Ok(output)
            }
        }

        deserializer.deserialize_map(ReadFieldVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::AuxArray;

    #[test]
    fn simple_aux_accepts_scalar_values() -> Result<(), Error> {
        assert_eq!(
            SimpleAux::try_from(Aux::Char(b'A'))?,
            SimpleAux::String("A".to_owned())
        );
        assert_eq!(SimpleAux::try_from(Aux::I8(-5))?, SimpleAux::I64(-5));
        assert_eq!(SimpleAux::try_from(Aux::U8(5))?, SimpleAux::I64(5));
        assert_eq!(SimpleAux::try_from(Aux::I16(-50))?, SimpleAux::I64(-50));
        assert_eq!(SimpleAux::try_from(Aux::U16(50))?, SimpleAux::I64(50));
        assert_eq!(SimpleAux::try_from(Aux::I32(-500))?, SimpleAux::I64(-500));
        assert_eq!(SimpleAux::try_from(Aux::U32(500))?, SimpleAux::I64(500));
        assert_eq!(
            SimpleAux::try_from(Aux::Float(1.5))?,
            SimpleAux::Double(1.5)
        );
        assert_eq!(
            SimpleAux::try_from(Aux::Double(2.5))?,
            SimpleAux::Double(2.5)
        );
        assert_eq!(
            SimpleAux::try_from(Aux::String("short"))?,
            SimpleAux::String("short".to_owned())
        );
        Ok(())
    }

    #[test]
    fn simple_aux_rejects_hex_byte_arrays() {
        let result = SimpleAux::try_from(Aux::HexByteArray("CAFE"));
        assert!(matches!(result, Err(Error::NotImplemented(_))));
    }

    #[test]
    fn simple_aux_from_impls() {
        assert_eq!(SimpleAux::from(-5i64), SimpleAux::I64(-5));
        assert_eq!(SimpleAux::from(1.25f64), SimpleAux::Double(1.25));
        assert_eq!(
            SimpleAux::from(String::from("hello")),
            SimpleAux::String(String::from("hello"))
        );
        assert_eq!(
            SimpleAux::from("world"),
            SimpleAux::String(String::from("world"))
        );
    }

    #[test]
    fn simple_aux_display() {
        assert_eq!(SimpleAux::I64(-5).to_string(), "-5");
        assert_eq!(SimpleAux::Double(1.25).to_string(), "1.25");
        assert_eq!(
            SimpleAux::String(String::from("hello")).to_string(),
            "hello"
        );
    }

    #[test]
    fn simple_aux_rejects_array_values() {
        let value = [1u8, 2u8, 3u8];
        let array: AuxArray<'_, u8> = (&value).into();
        let result = SimpleAux::try_from(Aux::ArrayU8(array));
        assert!(matches!(result, Err(Error::NotImplemented(_))));
    }

    #[test]
    fn read_field_new_and_wraps_btree_map() {
        let mut read_field = ReadField::new();
        assert_eq!(read_field.insert(*b"XU", 7i64), None);

        let mut expected = ReadField::default();
        assert_eq!(expected.insert(*b"XU", 7i64), None);

        assert_eq!(read_field, expected);
    }

    #[test]
    fn read_field_extend_owned_and_index() {
        let mut read_field = ReadField::default();
        read_field.extend([(*b"ts", 1i64), (*b"ul", 2i64)]);

        assert_eq!(read_field.get(&b"ts"[..]), Some(&SimpleAux::I64(1)));
        assert_eq!(read_field.get(&b"ul"[..]), Some(&SimpleAux::I64(2)));
    }

    #[test]
    fn read_field_clear_len_and_is_empty() {
        let mut read_field = ReadField::default();
        assert!(read_field.is_empty());
        assert_eq!(read_field.len(), 0);

        assert_eq!(read_field.insert(*b"ts", 1i64), None);
        assert_eq!(read_field.insert(*b"ul", 2i64), None);
        assert!(!read_field.is_empty());
        assert_eq!(read_field.len(), 2);

        read_field.clear();
        assert!(read_field.is_empty());
        assert_eq!(read_field.len(), 0);
    }

    #[test]
    fn read_field_get_contains_key_get_key_value_get_mut_and_entry() {
        let mut read_field = ReadField::new();

        for tag in [*b"ts", *b"ul", *b"ts"] {
            let _value = read_field
                .entry(tag)
                .and_modify(|value| {
                    if let &mut SimpleAux::I64(ref mut curr) = value {
                        *curr += 1;
                    }
                })
                .or_insert(SimpleAux::I64(1));
        }

        assert!(read_field.contains_key(&b"ts"[..]));
        assert!(read_field.contains_key(&b"ul"[..]));
        assert!(!read_field.contains_key(&b"db"[..]));
        assert_eq!(read_field.get(&b"ts"[..]), Some(&SimpleAux::I64(2)));
        assert_eq!(read_field.get(&b"db"[..]), None);
        assert_eq!(
            read_field.get_key_value(&b"ul"[..]),
            Some((b"ul", &SimpleAux::I64(1)))
        );

        let value = read_field
            .get_mut(&b"ul"[..])
            .expect("`ul` tag should exist in test setup");
        *value = SimpleAux::String(String::from("updated"));

        assert_eq!(read_field.get(&b"ts"[..]), Some(&SimpleAux::I64(2)));
        assert_eq!(
            read_field.get(&b"ul"[..]),
            Some(&SimpleAux::String(String::from("updated")))
        );
    }

    #[test]
    fn read_field_iter_keys_values_and_mut_variants() {
        let mut read_field = ReadField::new();
        assert_eq!(read_field.insert(*b"ts", String::from("hello")), None);
        assert_eq!(read_field.insert(*b"ul", String::from("goodbye")), None);

        let entries_before: Vec<([u8; 2], SimpleAux)> = read_field
            .iter()
            .map(|(key, value)| (*key, value.clone()))
            .collect();
        assert_eq!(
            entries_before,
            [
                (*b"ts", SimpleAux::String(String::from("hello"))),
                (*b"ul", SimpleAux::String(String::from("goodbye"))),
            ]
        );

        let keys: Vec<[u8; 2]> = read_field.keys().copied().collect();
        assert_eq!(keys, [*b"ts", *b"ul"]);

        let values_before: Vec<SimpleAux> = read_field.values().cloned().collect();
        assert_eq!(
            values_before,
            [
                SimpleAux::String(String::from("hello")),
                SimpleAux::String(String::from("goodbye")),
            ]
        );

        read_field.iter_mut().for_each(|(key, value)| {
            if key == b"ts" {
                *value = SimpleAux::String(String::from("hello?"));
            }
        });

        for value in read_field.values_mut() {
            if let &mut SimpleAux::String(ref mut text) = value {
                text.push('!');
            }
        }

        let values_after: Vec<SimpleAux> = read_field.values().cloned().collect();
        assert_eq!(
            values_after,
            [
                SimpleAux::String(String::from("hello?!")),
                SimpleAux::String(String::from("goodbye!")),
            ]
        );
    }

    #[test]
    fn read_field_extend_by_ref() {
        let entries = [
            (*b"ts", SimpleAux::I64(1)),
            (*b"ul", SimpleAux::String("blah".to_owned())),
        ];
        let mut read_field = ReadField::default();
        read_field.extend(entries.iter().map(|entry| (&entry.0, &entry.1)));

        assert_eq!(read_field.get(&b"ts"[..]), Some(&SimpleAux::I64(1)));
        assert_eq!(
            read_field.get(&b"ul"[..]),
            Some(&SimpleAux::String("blah".to_owned()))
        );
    }

    #[test]
    fn read_field_from_array_and_from_iter() {
        let from_array = ReadField::from([(*b"ts", 1i64), (*b"ul", 2i64)]);
        let from_iter = [(*b"ts", 1i64), (*b"ul", 2i64)]
            .into_iter()
            .collect::<ReadField>();

        assert_eq!(from_array, from_iter);
    }

    #[test]
    fn read_field_from_array_and_from_iter_strings() {
        let from_array =
            ReadField::from([(*b"ts", String::from("one")), (*b"ul", String::from("two"))]);
        let from_iter = [(*b"ts", "one"), (*b"ul", "two")]
            .into_iter()
            .collect::<ReadField>();

        assert_eq!(from_array, from_iter);
    }

    #[test]
    fn read_field_extend_owned_floats() {
        let mut read_field = ReadField::default();
        read_field.extend([(*b"ts", 1.5f64), (*b"ul", 2.5f64)]);

        assert_eq!(read_field.get(&b"ts"[..]), Some(&SimpleAux::Double(1.5)));
        assert_eq!(read_field.get(&b"ul"[..]), Some(&SimpleAux::Double(2.5)));
    }

    #[test]
    fn read_field_serializes_as_json_object() -> Result<(), Error> {
        let mut aux_map = BTreeMap::new();
        let _: Option<SimpleAux> = aux_map.insert(*b"ts", SimpleAux::I64(1));
        let _: Option<SimpleAux> = aux_map.insert(*b"ul", SimpleAux::String("blah".to_owned()));

        let json = serde_json::to_string(&ReadField(aux_map))?;
        assert_eq!(json, r#"{"ts":1,"ul":"blah"}"#);
        Ok(())
    }

    #[test]
    fn read_field_deserializes_from_json_object() -> Result<(), Error> {
        let read_field: ReadField = serde_json::from_str(r#"{"ts":1,"ul":"blah"}"#)?;

        let mut expected = BTreeMap::new();
        let _: Option<SimpleAux> = expected.insert(*b"ts", SimpleAux::I64(1));
        let _: Option<SimpleAux> = expected.insert(*b"ul", SimpleAux::String("blah".to_owned()));

        assert_eq!(read_field, ReadField(expected));
        Ok(())
    }

    #[test]
    fn read_field_serializes_all_three_datatypes() -> Result<(), Error> {
        let mut aux_map = BTreeMap::new();
        let _: Option<SimpleAux> = aux_map.insert(*b"ch", SimpleAux::String("A".to_owned()));
        let _: Option<SimpleAux> = aux_map.insert(*b"in", SimpleAux::I64(-42));
        let _: Option<SimpleAux> = aux_map.insert(*b"db", SimpleAux::Double(1.25));
        let _: Option<SimpleAux> = aux_map.insert(*b"st", SimpleAux::String("hello".to_owned()));

        let json = serde_json::to_string(&ReadField(aux_map))?;
        assert_eq!(json, r#"{"ch":"A","db":1.25,"in":-42,"st":"hello"}"#);
        Ok(())
    }

    #[test]
    fn read_field_deserializes_all_three_datatypes() -> Result<(), Error> {
        let read_field: ReadField =
            serde_json::from_str(r#"{"ch":"A","db":1.25,"in":-42,"st":"hello"}"#)?;

        let mut expected = BTreeMap::new();
        let _: Option<SimpleAux> = expected.insert(*b"ch", SimpleAux::String("A".to_owned()));
        let _: Option<SimpleAux> = expected.insert(*b"db", SimpleAux::Double(1.25));
        let _: Option<SimpleAux> = expected.insert(*b"in", SimpleAux::I64(-42));
        let _: Option<SimpleAux> = expected.insert(*b"st", SimpleAux::String("hello".to_owned()));

        assert_eq!(read_field, ReadField(expected));
        Ok(())
    }

    #[test]
    fn read_field_deserialization_rejects_non_two_byte_keys() {
        let result: Result<ReadField, serde_json::Error> = serde_json::from_str(r#"{"toolong":1}"#);
        let _: serde_json::Error = result.unwrap_err();
    }

    #[test]
    fn read_field_deserialization_rejects_duplicate_integer_keys() {
        let result: Result<ReadField, serde_json::Error> =
            serde_json::from_str(r#"{"ts":1,"ts":3}"#);
        let _: serde_json::Error = result.unwrap_err();
    }

    #[test]
    fn read_field_deserialization_rejects_duplicate_mixed_type_keys() {
        let result: Result<ReadField, serde_json::Error> =
            serde_json::from_str(r#"{"ts":1.11,"ts":"blah"}"#);
        let _: serde_json::Error = result.unwrap_err();
    }
}
