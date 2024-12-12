use std::{
    fs::File,
    io::{BufReader, BufWriter, Error, Read, Write},
    path::PathBuf,
    sync::Mutex,
};

use clap::Parser;
use log::{debug, info, LevelFilter};
use simplelog::{ColorChoice, CombinedLogger, TermLogger, TerminalMode};

static LOGGING_INITIALISED: Mutex<bool> = Mutex::new(false);

pub fn initialise_logging(log_level: LevelFilter) {
    let mut logging_initialised = LOGGING_INITIALISED.lock().unwrap();

    if !*logging_initialised {
        CombinedLogger::init(vec![TermLogger::new(
            log_level,
            Default::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto,
        )])
        .unwrap();

        info!("Logging initialised successfully");
        *logging_initialised = true;
    }
}

/// Upper case all genome characters and remove all non-ACGT characters.
#[derive(Parser, Debug)]
pub struct Config {
    /// The desired log level.
    #[clap(short, long, default_value = "Info")]
    log_level: LevelFilter,

    /// The input fasta file.
    #[clap(index = 1)]
    input: PathBuf,

    /// The output fasta file. Will be overwritten if it exists.
    #[clap(index = 2)]
    output: PathBuf,
}

enum FastaState {
    Init,
    RecordHeader {
        width: Option<usize>,
    },
    RecordHeaderLineBreak {
        width: Option<usize>,
    },
    RecordSequenceWithoutWidth {
        current_width: usize,
        current_output_row_width: usize,
    },
    RecordSequenceLineBreak {
        width: usize,
        current_output_row_width: usize,
    },
    RecordSequence {
        width: usize,
        current_output_row_width: usize,
    },
}

enum ReadResult<T> {
    Ok(T),
    Eof,
    Error(Error),
}

enum ReadOk<T> {
    Ok(T),
    Eof,
}

fn main() {
    let config = Config::parse();
    initialise_logging(config.log_level);
    debug!("{config:?}");

    info!("Opening input file: {:?}", config.input);
    let mut input = BufReader::new(File::open(&config.input).unwrap());
    info!("Opening output file: {:?}", config.output);
    let mut output = BufWriter::new(File::create(&config.output).unwrap());

    info!("Cleaning...");
    clean_fasta_file(&mut input, &mut output);

    // Manually calling drop here to ensure that "Done." is only printed after the files flushed and closed.
    drop(input);
    drop(output);
    info!("Done.");
}

fn clean_fasta_file(mut input: impl Read, mut output: impl Write) {
    let mut state = FastaState::Init;
    let mut buffer = Vec::new();

    loop {
        match state {
            FastaState::Init => match read_character(&mut input, &mut buffer).unwrap() {
                ReadOk::Ok(b'>') => {
                    state = FastaState::RecordHeader { width: None };
                    output.write_all(b">").unwrap();
                }
                ReadOk::Ok(character) => {
                    if !character.is_ascii_whitespace() {
                        panic!("Found non-whitespace character before first fasta record.");
                    }
                }
                ReadOk::Eof => break,
            },
            FastaState::RecordHeader { width } => {
                match read_character(&mut input, &mut buffer).unwrap() {
                    ReadOk::Ok(b'\n' | b'\r') => {
                        state = FastaState::RecordHeaderLineBreak { width };
                        output.write_all(b"\n").unwrap();
                    }
                    ReadOk::Ok(character) => output.write_all(&[character]).unwrap(),
                    ReadOk::Eof => break,
                }
            }
            FastaState::RecordHeaderLineBreak { width } => {
                match read_character(&mut input, &mut buffer).unwrap() {
                    ReadOk::Ok(b'\n' | b'\r') => { /* Ignore further line breaks */ }
                    ReadOk::Ok(b'>') => {
                        state = FastaState::RecordHeader { width };
                        output.write_all(b">").unwrap();
                    }
                    ReadOk::Ok(character) => {
                        let character = character.to_ascii_uppercase();
                        if let Some(width) = width {
                            state = FastaState::RecordSequence {
                                width,
                                current_output_row_width: 1,
                            };
                        } else {
                            state = FastaState::RecordSequenceWithoutWidth {
                                current_width: 1,
                                current_output_row_width: 1,
                            };
                        }
                        output.write_all(&[character]).unwrap();
                    }
                    ReadOk::Eof => break,
                }
            }
            FastaState::RecordSequenceWithoutWidth {
                mut current_width,
                mut current_output_row_width,
            } => match read_character(&mut input, &mut buffer).unwrap() {
                ReadOk::Ok(b'\n' | b'\r') => {
                    debug_assert!(current_width > 0);
                    if current_output_row_width == current_width {
                        output.write_all(b"\n").unwrap();
                        current_output_row_width = 0;
                    }
                    debug!("Found fasta line width {current_width}");
                    state = FastaState::RecordSequenceLineBreak {
                        width: current_width,
                        current_output_row_width,
                    };
                }
                ReadOk::Ok(b'>') => panic!("Encountered '>' within sequence."),
                ReadOk::Ok(character) => {
                    let character = character.to_ascii_uppercase();
                    current_width += 1;

                    if matches!(character, b'A' | b'C' | b'G' | b'T') {
                        current_output_row_width += 1;
                        output.write_all(&[character]).unwrap();
                    }

                    state = FastaState::RecordSequenceWithoutWidth {
                        current_width,
                        current_output_row_width,
                    };
                }
                ReadOk::Eof => {
                    output.write_all(b"\n").unwrap();
                    break;
                }
            },
            FastaState::RecordSequenceLineBreak {
                width,
                mut current_output_row_width,
            } => match read_character(&mut input, &mut buffer).unwrap() {
                ReadOk::Ok(b'\n' | b'\r') => { /* Ignore further line breaks */ }
                ReadOk::Ok(b'>') => {
                    state = FastaState::RecordHeader { width: Some(width) };
                    debug_assert!(current_output_row_width <= width);
                    if current_output_row_width > 0 {
                        output.write_all(b"\n").unwrap();
                    }
                    output.write_all(b">").unwrap();
                }
                ReadOk::Ok(character) => {
                    let character = character.to_ascii_uppercase();

                    debug_assert!(current_output_row_width <= width);
                    if current_output_row_width == width {
                        output.write_all(b"\n").unwrap();
                        current_output_row_width = 0;
                    }

                    if matches!(character, b'A' | b'C' | b'G' | b'T') {
                        current_output_row_width += 1;
                        let character = character.to_ascii_uppercase();
                        output.write_all(&[character]).unwrap();
                    }

                    state = FastaState::RecordSequence {
                        width,
                        current_output_row_width,
                    };
                }
                ReadOk::Eof => {
                    output.write_all(b"\n").unwrap();
                    break;
                }
            },
            FastaState::RecordSequence {
                width,
                mut current_output_row_width,
            } => match read_character(&mut input, &mut buffer).unwrap() {
                ReadOk::Ok(b'\n' | b'\r') => {
                    debug_assert!(current_output_row_width <= width);
                    if current_output_row_width == width {
                        output.write_all(b"\n").unwrap();
                        current_output_row_width = 0;
                    }
                    state = FastaState::RecordSequenceLineBreak {
                        width,
                        current_output_row_width,
                    };
                }
                ReadOk::Ok(b'>') => panic!("Encountered '>' within sequence."),
                ReadOk::Ok(character) => {
                    let character = character.to_ascii_uppercase();

                    debug_assert!(current_output_row_width <= width);
                    if current_output_row_width == width {
                        output.write_all(b"\n").unwrap();
                        current_output_row_width = 0;
                    }

                    if matches!(character, b'A' | b'C' | b'G' | b'T') {
                        current_output_row_width += 1;
                        let character = character.to_ascii_uppercase();
                        output.write_all(&[character]).unwrap();
                    }

                    state = FastaState::RecordSequence {
                        width,
                        current_output_row_width,
                    };
                }
                ReadOk::Eof => {
                    output.write_all(b"\n").unwrap();
                    break;
                }
            },
        }
    }
}

fn read_buffer(reader: &mut impl Read, buffer: &mut Vec<u8>, length: usize) -> ReadResult<()> {
    buffer.resize(length, 0);
    match reader.read_exact(buffer) {
        Ok(()) => ReadResult::Ok(()),
        Err(error) => match error.kind() {
            std::io::ErrorKind::UnexpectedEof => ReadResult::Eof,
            _ => ReadResult::Error(error),
        },
    }
}

fn read_character(reader: &mut impl Read, buffer: &mut Vec<u8>) -> ReadResult<u8> {
    match read_buffer(reader, buffer, 1) {
        ReadResult::Ok(()) => ReadResult::Ok(buffer[0]),
        ReadResult::Eof => ReadResult::Eof,
        ReadResult::Error(error) => ReadResult::Error(error),
    }
}

impl<T> ReadResult<T> {
    pub fn unwrap(self) -> ReadOk<T> {
        match self {
            ReadResult::Ok(value) => ReadOk::Ok(value),
            ReadResult::Eof => ReadOk::Eof,
            ReadResult::Error(error) => panic!("read error: {error}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use log::{debug, LevelFilter};

    use crate::{clean_fasta_file, initialise_logging};

    #[test]
    fn test() {
        initialise_logging(LevelFilter::Debug);
        test_file(
            b"\r>WGCaC\n\nAACCcxXAA\naacc\n.ef34\nCGG\ntgtcgcgtagcgtgatcgtgtagtcgtag\r.\r>f\nTTT",
            b">WGCaC\nAACCCAAAA\nCCCGGTGTC\nGCGTAGCGT\nGATCGTGTA\nGTCGTAG\n>f\nTTT\n",
        );
    }

    fn test_file(input: &[u8], expected_output: &[u8]) {
        debug!("input:\n{}", String::from_utf8_lossy(input));
        let mut output = Vec::new();

        clean_fasta_file(input, &mut output);
        assert_eq!(
            output,
            expected_output,
            "actual:\n{}\nexpected:\n{}",
            String::from_utf8_lossy(&output),
            String::from_utf8_lossy(expected_output),
        );
    }
}
