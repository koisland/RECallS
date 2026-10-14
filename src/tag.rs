use std::{
    collections::{HashMap, HashSet},
    fs::File,
    num::NonZero,
    path::Path,
};

use noodles::{
    bam, bgzf,
    core::{Position, Region},
    sam::alignment::{
        RecordBuf, io::Write, record::data::field::Tag, record_buf::data::field::Value,
    },
};

use crate::events::{SV, SVType};

pub fn create_tagged_bam(bam: &Path, out_bam: &Path, events: &HashMap<String, Vec<SV>>) -> eyre::Result<()> {
    let mut fh = bam::io::indexed_reader::Builder::default().build_from_path(bam)?;
    let header = fh.read_header()?;

    let dst_file = File::create(out_bam)?;
    let encoder =
        bgzf::io::MultithreadedWriter::with_worker_count(NonZero::new(4).unwrap(), dst_file);
    let mut fh_out = bam::io::Writer::from(encoder);
    fh_out.write_header(&header)?;

    // Partition reads by event
    let (inv_events, del_events) = events.values().fold(
        (HashSet::new(), HashSet::new()),
        |(mut inv_events, mut del_events), svs| {
            for (sv_type, reads) in svs.iter().map(|sv| (sv.typ.clone(), sv.reads())) {
                match sv_type {
                    SVType::Deletion => {
                        del_events.extend(reads);
                    }
                    SVType::Inversion => {
                        inv_events.extend(reads);
                    }
                }
            }
            (inv_events, del_events)
        },
    );

    // TODO: Parallelize
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
            let tag_value = if inv_events.contains(rname) {
                Some("Inversion")
            } else if del_events.contains(rname) {
                Some("Deletion")
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

    // Create bam index
    std::mem::drop(fh_out);
    let index = bam::fs::index(bam)?;
    let mut bai = bam.to_path_buf();
    bai.add_extension("bai");
    bam::bai::fs::write(bai, &index)?;

    Ok(())
}
