use itertools::Itertools;
use rust_lapper::Interval;

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

pub fn merge_intervals<I, T>(
    intervals: I,
    dst: usize,
    fn_cmp: impl Fn(&Interval<usize, T>, &Interval<usize, T>) -> bool,
    fn_reducer: impl Fn(&Interval<usize, T>, &Interval<usize, T>) -> T,
) -> Vec<Interval<usize, T>>
where
    I: Iterator<Item = Interval<usize, T>>,
    T: Clone + Eq + Send + Sync,
{
    let mut iter_intervals = intervals.sorted_by(|a, b| a.start.cmp(&b.start));
    let Some(itv_first) = iter_intervals.next() else {
        return vec![];
    };
    let mut merged = vec![itv_first];

    for itv in iter_intervals {
        if let Some(prev) = merged.last_mut().filter(|prev| {
            let dst_between = itv.start.saturating_sub(prev.stop);
            let added_check = fn_cmp(prev, &itv);
            (dst_between <= dst) & added_check
        }) {
            prev.stop = itv.stop;
            prev.val = fn_reducer(prev, &itv)
        } else {
            merged.push(itv);
        }
    }
    merged
}
