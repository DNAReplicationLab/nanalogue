#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Public interval intersection coverage, including empty and reversed ranges.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::Intersects as _;
    use std::collections::BTreeSet;
    use std::ops::Range;

    #[test]
    fn empty_ranges_inside_nonempty_ranges_do_not_intersect() {
        // An interior empty range passes both endpoint comparisons. It must
        // still be rejected by the explicit emptiness check, on either side.
        let enclosing = 3u32..11;
        for point in [3, 4, 7, 10, 11] {
            let empty = point..point;
            assert!(
                !empty.intersects(&enclosing),
                "empty receiver at {point} cannot intersect a populated range"
            );
            assert!(
                !enclosing.intersects(&empty),
                "empty argument at {point} cannot intersect a populated range"
            );
            assert!(
                !empty.intersects(&empty),
                "an empty range cannot intersect itself"
            );
        }
    }

    #[test]
    fn reversed_ranges_are_empty_even_when_endpoints_overlap() {
        // These values are constructible through the public Range type, so
        // unlike an invalid OrdPair they are reachable without unsafe code.
        let enclosing = 2u32..13;
        for (start, end) in [(9, 5), (12, 3), (13, 2), (u32::MAX, 0)] {
            let reversed = Range { start, end };
            assert!(reversed.is_empty(), "reversed endpoints describe no values");
            assert!(
                !reversed.intersects(&enclosing),
                "reversed receiver {reversed:?} must be rejected"
            );
            assert!(
                !enclosing.intersects(&reversed),
                "reversed argument {reversed:?} must be rejected"
            );
        }
    }

    #[test]
    fn respects_half_open_endpoints_without_coordinate_arithmetic() {
        // Asymmetric overlaps distinguish intersection from containment.
        // Near MAX, converting an exclusive endpoint to an inclusive one or
        // adding one to a coordinate is unnecessary and risks overflow.
        for (left, right, expected) in [
            (0u32..1, 1..4, false),
            (3..7, 8..12, false),
            (3..8, 7..12, true),
            (3..12, 5..8, true),
            (5..8, 5..8, true),
            (0..u32::MAX, 0..1, true),
            (0..u32::MAX, 0xFFFF_FFFE..u32::MAX, true),
            (0xFFFF_FFFD..u32::MAX, 0xFFFF_FFFE..u32::MAX, true),
            (0..u32::MAX, u32::MAX..u32::MAX, false),
            (0..0, 0..u32::MAX, false),
        ] {
            assert_eq!(
                left.intersects(&right),
                expected,
                "incorrect intersection for {left:?} and {right:?}"
            );
            assert_eq!(
                right.intersects(&left),
                expected,
                "intersection must be symmetric for {left:?} and {right:?}"
            );
        }
    }

    #[test]
    fn small_domain_matches_independent_set_intersection() {
        // Enumerating actual members gives an oracle independent of the
        // production endpoint comparisons and explicit emptiness guards.
        // Every ordering is included, not just well-formed intervals.
        let ranges: Vec<Range<u32>> = (0..=5)
            .flat_map(|start| (0..=5).map(move |end| start..end))
            .collect();
        for left in &ranges {
            let left_members: BTreeSet<u32> = left.clone().collect();
            for right in &ranges {
                let right_members: BTreeSet<u32> = right.clone().collect();
                let expected = !left_members.is_disjoint(&right_members);
                assert_eq!(
                    left.intersects(right),
                    expected,
                    "endpoint implementation differs from set intersection: {left:?}, {right:?}"
                );
            }
        }
    }
}
