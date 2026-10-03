use std::{collections::HashMap, fs::File, io::BufWriter, num::NonZero, path::Path};

use noodles::{
    bam, bgzf, core::{Position, Region}, sam::alignment::{
        RecordBuf, io::Write, record::data::field::Tag, record_buf::data::field::Value,
    },
};
// use rust_lapper::Interval;

use crate::events::{MismatchSignal, SupplSignal};

pub fn tag_bam(
    bam: &Path,
    out_bam: &Path,
    inv_events: &HashMap<String, Vec<MismatchSignal>>,
    del_events: &HashMap<String, Vec<SupplSignal>>,
) -> eyre::Result<()> {
    let mut fh = bam::io::indexed_reader::Builder::default().build_from_path(bam)?;
    let header = fh.read_header()?;

    let dst_file = File::create(out_bam)?;
    let encoder = bgzf::io::MultithreadedWriter::with_worker_count(NonZero::new(4).unwrap(), dst_file);
    let mut fh_out = bam::io::Writer::from(encoder);
    let mut new_header = header.clone();
    // new_header.add_comment("@PG     ID:RECallS     PN:RECallS     VN:0.0.1");
    fh_out.write_header(&new_header)?;

    for region in header
        .reference_sequences()
        .into_iter()
        .map(move |(ctg, ref_seq)| {
            let length: usize = ref_seq.length().get();
            Region::new(
                ctg.to_owned(),
                Position::new(1).unwrap()..=Position::new(length).unwrap(),
            )
        })
    {
        let query = fh.query(&header, &region)?;
        for rec in query.records().flatten() {
            let rname = str::from_utf8(rec.name().unwrap())?;
            // let (rst, rend) = (
            //     rec.alignment_start().unwrap().map(|p| p.get())?,
            //     rec.alignment_end().unwrap().map(|p| p.get())?,
            // );
            let tag_value = if inv_events.contains_key(rname) {
                Some("inv")
            } else if del_events.contains_key(rname) {
                Some("del")
            } else {
                None
            };

            match tag_value {
                Some(tag_value) => {
                    let mut rec_buf = RecordBuf::try_from_alignment_record(&header, &rec)?;
                    let tags = rec_buf.data_mut();
                    tags.insert(Tag::new(b'E', b'V'), Value::String(tag_value.into()));
                    fh_out.write_alignment_record(&header, &rec_buf)?;
                }
                None => {
                    // fh_out.write_record(&header, &rec)?;
                }
            }
        }
    }
    Ok(())
}
