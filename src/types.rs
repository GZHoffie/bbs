// This file contains multiple implementations from sylph (https://github.com/bluenote-1577/sylph). Below is their license.

/*
MIT License

Copyright (c) 2023 Jim Shaw

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

use std::collections::HashMap;

pub type Kmer = u64;

// use 3 bits per nucleotide to allow for a special end marker (0b100).
pub const NUM_BITS_PER_NUC: usize = 3;
pub const END_MARKER: u8 = 0b100;

/**
 * A lookup table to convert a byte to a 2-bit sequence.
 * Adopted from https://github.com/bluenote-1577/sylph/blob/main/src/types.rs
 * 
 * A/a -> 0
 * C/c -> 1
 * G/g -> 2
 * T/t -> 3, U/u -> 3
 */
pub const BYTE_TO_SEQ: [u8; 256] = {
    let mut arr = [0u8; 256];

    arr[b'A' as usize] = 0;
    arr[b'C' as usize] = 1;
    arr[b'G' as usize] = 2;
    arr[b'T' as usize] = 3;
    arr[b'U' as usize] = 3;

    arr[b'a' as usize] = 0;
    arr[b'c' as usize] = 1;
    arr[b'g' as usize] = 2;
    arr[b't' as usize] = 3;
    arr[b'u' as usize] = 3;

    arr[b'$' as usize] = END_MARKER; // special end marker, 0b100.

    // [TODO] Handle other non-canonical bases if needed.

    arr
};

pub const SEQ_TO_BYTE: [u8; 5] = [b'A', b'C', b'G', b'T', b'$'];