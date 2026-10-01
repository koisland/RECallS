#[inline]
// Adapted from https://github.com/chaimleib/intervaltree/blob/1bc406e1f441577c4e421fc51aba2ab67fbd97fb/intervaltree/interval.py#L56
pub fn overlap_length(a_first: usize, a_last: usize, b_first: usize, b_last: usize) -> usize {
    // No overlap
    if !(a_first < b_last && a_last > b_first) {
        return 0;
    }
    // a  |---|
    // b |---|
    let max_st = std::cmp::max(a_first, b_first);
    let min_end = std::cmp::min(a_last, b_last);
    min_end - max_st
}
