//! FilterByRefCoords trait for filtering by coordinates on the reference genome
//! Provides interface for coordinate-based filtering operations

/// Implements filter by coordinates on the reference genome.
pub trait FilterByRefCoords {
    /// filters by reference position i.e. all pos such that start <= pos < end
    /// are retained. does not use contig in filtering.
    fn filter_by_ref_pos(&mut self, _: i64, _: i64) {
        todo!()
    }
}
