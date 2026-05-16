//! DNA complement and reverse-complement helpers.
//!
//! This module includes code adapted from the `bio` crate's DNA alphabet
//! utilities. See `THIRD_PARTY_NOTICES.md` for the upstream license text.

use std::array;
use std::borrow::Borrow;
use std::sync::LazyLock;

/// Complement lookup table adapted from the `bio` crate's DNA alphabet helpers.
static COMPLEMENT: LazyLock<[u8; 256]> = LazyLock::new(|| {
    let mut comp = array::from_fn(|v| u8::try_from(v).expect("array indices 0..=255 fit in u8"));
    b"AGCTYRWSKMDVHBN"
        .iter()
        .zip(b"TCGARYWSMKHBDVN".iter())
        .for_each(|(&a, &b)| {
            *comp
                .get_mut(a as usize)
                .expect("ASCII nucleotide bytes are valid lookup indices") = b;
            *comp
                .get_mut(a as usize + 32)
                .expect("ASCII lowercase nucleotide bytes are valid lookup indices") = b + 32;
        });
    comp
});

/// Return the DNA complement of a byte while preserving `bio`'s permissive behavior.
///
/// Unknown bytes are returned unchanged. IUPAC ambiguity codes and lowercase
/// bases are supported.
///
/// # Panics
///
/// Panics only if a `u8` fails to index the 256-byte lookup table, which is
/// impossible.
#[inline]
#[must_use]
pub fn complement(a: u8) -> u8 {
    COMPLEMENT
        .get(a as usize)
        .copied()
        .expect("u8 values always index into the 256-byte complement table")
}

/// Return the reverse complement of an input byte sequence.
///
/// The fast path (`.collect()`) is only taken when the iterator's reported
/// upper-bound length is within a safe pre-reservation budget. Larger
/// reported lengths are refused and return an empty `Vec` instead.
#[must_use]
pub fn revcomp<C, T>(text: T) -> Vec<u8>
where
    C: Borrow<u8>,
    T: IntoIterator<Item = C>,
    T::IntoIter: DoubleEndedIterator,
{
    /// Upper bound on the size-hint we are willing to pre-allocate (3 GiB).
    /// Anything reporting more is refused.
    const MAX_PRE_RESERVE: usize = 3usize * 1024 * 1024 * 1024;

    let iter = text.into_iter().rev().map(|a| complement(*a.borrow()));
    let (lo, hi) = iter.size_hint();
    let upper = hi.unwrap_or(lo);

    if upper <= MAX_PRE_RESERVE {
        iter.collect()
    } else {
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_preserves_iupac_behavior() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'c'), b'g');
        assert_eq!(complement(b'N'), b'N');
        assert_eq!(complement(b'Y'), b'R');
        assert_eq!(complement(b's'), b's');
        assert_eq!(complement(b'='), b'=');
        assert_eq!(complement(b'X'), b'X');
    }

    #[test]
    fn complement_table_has_all_u8_entries() {
        let table = &*COMPLEMENT;
        assert_eq!(table.len(), usize::from(u8::MAX) + 1);

        for byte in u8::MIN..=u8::MAX {
            assert!(table.get(usize::from(byte)).is_some());
        }
    }

    #[test]
    fn revcomp_preserves_iupac_behavior() {
        assert_eq!(revcomp(b"ACGTN"), b"NACGT");
        assert_eq!(revcomp(b"GaTtaCA"), b"TGtaAtC");
        assert_eq!(revcomp(b"AGCTYRWSKMDVHBN"), b"NVDBHKMSWYRAGCT");
    }

    /// Iterators with `size_hint` exceeding the pre-reserve budget return an
    /// empty `Vec`.
    #[test]
    fn revcomp_does_not_abort_on_adversarial_size_hint() {
        let it = std::iter::repeat_n(0u8, usize::MAX);
        let out = revcomp(it);
        assert!(
            out.is_empty(),
            "adversarial size_hint must produce an empty Vec, got {}",
            out.len(),
        );
    }
}
