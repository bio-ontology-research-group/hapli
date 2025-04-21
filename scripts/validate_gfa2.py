#!/usr/bin/env python3

import argparse
import sys
import re
from typing import List, Set, Tuple

# --- GFA2 Regex Definitions ---
# Basic types
ID_RE = re.compile(r'^[!-~]+$')
OPT_ID_RE = re.compile(r'^(\*|[!-~]+)$')
INT_RE = re.compile(r'^-?[0-9]+$')
UINT_RE = re.compile(r'^[0-9]+$') # Unsigned integer for lengths, counts
POS_RE = re.compile(r'^-?[0-9]+\$?$') # Position (integer optionally ending with $)
REF_RE = re.compile(r'^[!-~]+[+-]$') # Reference (ID followed by + or -)
SEQUENCE_RE = re.compile(r'^(\*|[A-Za-z=.]+)$') # GFA2 Spec says [!-~]+, but often restricted like SAM
# Allowing full [!-~]+ for sequence as per strict spec interpretation:
# SEQUENCE_RE = re.compile(r'^(\*|[!-~]+)$')
TAG_RE = re.compile(r'^[A-Za-z0-9]{2}:[ABHJZif]:[ -~]*$')

# Alignment types
CIGAR_RE = re.compile(r'^(\*|([0-9]+[MDIPN=X])+)?$') # Added N=X based on common usage, * is valid empty alignment
TRACE_RE = re.compile(r'^-?[0-9]+(,-?[0-9]+)*$')
# Alignment allows *, trace, or CIGAR. '*' is handled separately.
ALIGNMENT_CONTENT_RE = re.compile(f'^({TRACE_RE.pattern}|{CIGAR_RE.pattern})$')

# --- Helper Functions ---

def _validate_tag(tag: str, line_num: int, record_type: str, field_idx: int, errors: List[str]):
    """Validates a single optional tag."""
    if not TAG_RE.match(tag):
        errors.append(f"L{line_num} ({record_type}): Invalid tag format '{tag}' at field index {field_idx}.")
        return False
    # Basic format is valid, further type/value checks could be added here if needed
    # e.g., check if value matches type 'i' for integer, 'f' for float etc.
    return True

def _validate_id(id_val: str, line_num: int, record_type: str, field_idx: int, errors: List[str], id_type: str = "ID") -> bool:
    """Validates an ID field."""
    if not ID_RE.match(id_val):
        errors.append(f"L{line_num} ({record_type}): Invalid {id_type} '{id_val}' at field index {field_idx}. Must match [!-~]+.")
        return False
    return True

def _validate_opt_id(id_val: str, line_num: int, record_type: str, field_idx: int, errors: List[str], id_type: str = "ID") -> bool:
    """Validates an optional ID field (* or ID)."""
    if not OPT_ID_RE.match(id_val):
        errors.append(f"L{line_num} ({record_type}): Invalid optional {id_type} '{id_val}' at field index {field_idx}. Must be '*' or match [!-~]+.")
        return False
    return True

def _validate_ref(ref_val: str, line_num: int, record_type: str, field_idx: int, errors: List[str], ref_type: str = "reference") -> bool:
    """Validates a reference field (ID+ or ID-)."""
    if not REF_RE.match(ref_val):
        errors.append(f"L{line_num} ({record_type}): Invalid {ref_type} '{ref_val}' at field index {field_idx}. Must be ID followed by '+' or '-'.")
        return False
    # Could add check here if the base ID exists in defined_ids if strict checking is needed
    return True

def _validate_pos(pos_val: str, line_num: int, record_type: str, field_idx: int, errors: List[str], pos_type: str = "position") -> bool:
    """Validates a position field (integer optionally ending in $)."""
    if not POS_RE.match(pos_val):
        errors.append(f"L{line_num} ({record_type}): Invalid {pos_type} '{pos_val}' at field index {field_idx}. Must be integer, optionally ending with '$'.")
        return False
    # Check if value is reasonable (e.g. >= 0 if $ is not present?) - depends on context
    val = int(pos_val.rstrip('$'))
    if val < 0:
         errors.append(f"L{line_num} ({record_type}): Position '{pos_val}' at field index {field_idx} cannot be negative.")
         return False
    return True

def _validate_int(int_val: str, line_num: int, record_type: str, field_idx: int, errors: List[str], int_type: str = "integer") -> bool:
    """Validates an integer field."""
    if not INT_RE.match(int_val):
        errors.append(f"L{line_num} ({record_type}): Invalid {int_type} '{int_val}' at field index {field_idx}. Must be an integer.")
        return False
    return True

def _validate_uint(uint_val: str, line_num: int, record_type: str, field_idx: int, errors: List[str], uint_type: str = "integer") -> bool:
    """Validates a non-negative integer field."""
    if not UINT_RE.match(uint_val):
        errors.append(f"L{line_num} ({record_type}): Invalid {uint_type} '{uint_val}' at field index {field_idx}. Must be a non-negative integer.")
        return False
    return True

def _validate_alignment(aln_val: str, line_num: int, record_type: str, field_idx: int, errors: List[str]) -> bool:
    """Validates an alignment field (*, CIGAR, or trace)."""
    if aln_val == '*':
        return True
    # Check CIGAR first as it's more specific
    if CIGAR_RE.match(aln_val):
        # Further check if it's *not* just '*' which is already handled
        if aln_val != '*':
             # Could add CIGAR semantic validation here (e.g., non-zero lengths)
             return True
    # Check trace
    if TRACE_RE.match(aln_val):
        return True

    errors.append(f"L{line_num} ({record_type}): Invalid alignment format '{aln_val}' at field index {field_idx}. Must be '*', valid CIGAR, or trace.")
    return False

# --- Record Type Validation Functions ---

def validate_header(fields: List[str], line_num: int, errors: List[str]):
    # H {VN:Z:2.0} {TS:i:<trace spacing>} <tag>*
    # Spec implies VN and TS are optional tags, not fixed fields.
    # We only check the record type 'H'. All other fields must be valid tags.
    if len(fields) < 1:
        errors.append(f"L{line_num} (H): Header line is empty.") # Should not happen if line stripping works
        return

    has_vn = False
    for i, tag in enumerate(fields[1:]):
        if not _validate_tag(tag, line_num, 'H', i + 1, errors):
            continue # Skip further checks on this invalid tag
        tag_parts = tag.split(':')
        tag_name = tag_parts[0]
        tag_type = tag_parts[1]
        tag_value = tag_parts[2]

        if tag_name == "VN":
            has_vn = True
            if tag_type != 'Z':
                 errors.append(f"L{line_num} (H): VN tag must have type 'Z', found '{tag_type}'.")
            if tag_value != "2.0":
                 errors.append(f"L{line_num} (H): VN:Z tag value must be '2.0', found '{tag_value}'.")
        elif tag_name == "TS":
            if tag_type != 'i':
                 errors.append(f"L{line_num} (H): TS tag must have type 'i', found '{tag_type}'.")
            if not INT_RE.match(tag_value): # Check if value is integer
                 errors.append(f"L{line_num} (H): TS:i tag value must be an integer, found '{tag_value}'.")
            # else: could check if TS value is positive?

    # GFA2 spec doesn't strictly require H line or VN tag, but it's good practice.
    # if not has_vn:
    #    errors.append(f"L{line_num} (H): Missing required VN:Z:2.0 tag.")


def validate_segment(fields: List[str], line_num: int, errors: List[str], defined_ids: Set[str]):
    # S <sid:id> <slen:int> <sequence> <tag>*
    min_fields = 4
    if len(fields) < min_fields:
        errors.append(f"L{line_num} (S): Expected at least {min_fields} fields, got {len(fields)}.")
        return

    # Field 1: sid
    sid = fields[1]
    if _validate_id(sid, line_num, 'S', 1, errors, id_type="segment ID"):
        if sid in defined_ids:
            errors.append(f"L{line_num} (S): Duplicate ID definition '{sid}'.")
        else:
            defined_ids.add(sid)

    # Field 2: slen
    slen_str = fields[2]
    slen = -1 # Sentinel for invalid length
    if _validate_uint(slen_str, line_num, 'S', 2, errors, uint_type="segment length"):
        slen = int(slen_str)

    # Field 3: sequence
    sequence = fields[3]
    if not SEQUENCE_RE.match(sequence):
         errors.append(f"L{line_num} (S): Invalid sequence characters found in '{sequence}' at field index 3. Must be '*' or match [A-Za-z=.]+ (or [!-~]+ depending on interpretation).")
    elif sequence != '*' and slen != -1 and len(sequence) != slen:
         errors.append(f"L{line_num} (S): Sequence length {len(sequence)} does not match declared length {slen}.")
    elif sequence == '*' and slen == -1:
         errors.append(f"L{line_num} (S): Sequence is '*' but segment length is invalid or missing.")
    # '*' sequence with valid slen >= 0 is okay.

    # Fields 4+: tags
    for i, tag in enumerate(fields[min_fields:]):
        _validate_tag(tag, line_num, 'S', i + min_fields, errors)


def validate_fragment(fields: List[str], line_num: int, errors: List[str], defined_ids: Set[str]):
    # F <sid:id> <external:ref> <sbeg:pos> <send:pos> <fbeg:pos> <fend:pos> <alignment> <tag>*
    min_fields = 8
    if len(fields) < min_fields:
        errors.append(f"L{line_num} (F): Expected at least {min_fields} fields, got {len(fields)}.")
        return

    # Field 1: sid (Segment ID this fragment belongs to)
    sid = fields[1]
    if not _validate_id(sid, line_num, 'F', 1, errors, id_type="segment ID"):
        pass # Error already added
    # elif sid not in defined_ids: # Check if segment ID exists - optional semantic check
    #     errors.append(f"L{line_num} (F): Referenced segment ID '{sid}' not defined.")

    # Field 2: external (Reference sequence name with orientation)
    _validate_ref(fields[2], line_num, 'F', 2, errors, ref_type="external reference")

    # Field 3: sbeg (Segment begin position)
    _validate_pos(fields[3], line_num, 'F', 3, errors, pos_type="segment begin")

    # Field 4: send (Segment end position)
    _validate_pos(fields[4], line_num, 'F', 4, errors, pos_type="segment end")

    # Field 5: fbeg (Fragment begin on external reference)
    _validate_pos(fields[5], line_num, 'F', 5, errors, pos_type="fragment begin")

    # Field 6: fend (Fragment end on external reference)
    _validate_pos(fields[6], line_num, 'F', 6, errors, pos_type="fragment end")

    # Field 7: alignment
    _validate_alignment(fields[7], line_num, 'F', 7, errors)

    # Fields 8+: tags
    for i, tag in enumerate(fields[min_fields:]):
        _validate_tag(tag, line_num, 'F', i + min_fields, errors)


def validate_edge(fields: List[str], line_num: int, errors: List[str], defined_ids: Set[str]):
    # E <eid:opt_id> <sid1:ref> <sid2:ref> <beg1:pos> <end1:pos> <beg2:pos> <end2:pos> <alignment> <tag>*
    min_fields = 9
    if len(fields) < min_fields:
        errors.append(f"L{line_num} (E): Expected at least {min_fields} fields, got {len(fields)}.")
        return

    # Field 1: eid (Optional Edge ID)
    eid = fields[1]
    if _validate_opt_id(eid, line_num, 'E', 1, errors, id_type="edge ID"):
        if eid != '*' and eid in defined_ids:
            errors.append(f"L{line_num} (E): Duplicate ID definition '{eid}'.")
        elif eid != '*':
            defined_ids.add(eid)

    # Field 2: sid1 (First segment reference)
    _validate_ref(fields[2], line_num, 'E', 2, errors, ref_type="segment 1 reference")
    # sid1_base = fields[2][:-1]
    # if sid1_base not in defined_ids: errors.append(...) # Optional semantic check

    # Field 3: sid2 (Second segment reference)
    _validate_ref(fields[3], line_num, 'E', 3, errors, ref_type="segment 2 reference")
    # sid2_base = fields[3][:-1]
    # if sid2_base not in defined_ids: errors.append(...) # Optional semantic check

    # Field 4: beg1 (Position on segment 1)
    _validate_pos(fields[4], line_num, 'E', 4, errors, pos_type="segment 1 begin")

    # Field 5: end1 (Position on segment 1)
    _validate_pos(fields[5], line_num, 'E', 5, errors, pos_type="segment 1 end")

    # Field 6: beg2 (Position on segment 2)
    _validate_pos(fields[6], line_num, 'E', 6, errors, pos_type="segment 2 begin")

    # Field 7: end2 (Position on segment 2)
    _validate_pos(fields[7], line_num, 'E', 7, errors, pos_type="segment 2 end")

    # Field 8: alignment
    _validate_alignment(fields[8], line_num, 'E', 8, errors)

    # Fields 9+: tags
    for i, tag in enumerate(fields[min_fields:]):
        _validate_tag(tag, line_num, 'E', i + min_fields, errors)


def validate_gap(fields: List[str], line_num: int, errors: List[str], defined_ids: Set[str]):
    # G <gid:opt_id> <sid1:ref> <sid2:ref> <dist:int> <var:int | *> <tag>*
    min_fields = 6
    if len(fields) < min_fields:
        errors.append(f"L{line_num} (G): Expected at least {min_fields} fields, got {len(fields)}.")
        return

    # Field 1: gid (Optional Gap ID)
    gid = fields[1]
    if _validate_opt_id(gid, line_num, 'G', 1, errors, id_type="gap ID"):
        if gid != '*' and gid in defined_ids:
            errors.append(f"L{line_num} (G): Duplicate ID definition '{gid}'.")
        elif gid != '*':
            defined_ids.add(gid)

    # Field 2: sid1 (First segment reference)
    _validate_ref(fields[2], line_num, 'G', 2, errors, ref_type="segment 1 reference")
    # sid1_base = fields[2][:-1]
    # if sid1_base not in defined_ids: errors.append(...) # Optional semantic check

    # Field 3: sid2 (Second segment reference)
    _validate_ref(fields[3], line_num, 'G', 3, errors, ref_type="segment 2 reference")
    # sid2_base = fields[3][:-1]
    # if sid2_base not in defined_ids: errors.append(...) # Optional semantic check

    # Field 4: dist (Distance)
    _validate_int(fields[4], line_num, 'G', 4, errors, int_type="distance")

    # Field 5: var (Variance)
    var_val = fields[5]
    if var_val != '*' and not INT_RE.match(var_val):
        errors.append(f"L{line_num} (G): Invalid variance '{var_val}' at field index 5. Must be '*' or integer.")

    # Fields 6+: tags
    for i, tag in enumerate(fields[min_fields:]):
        _validate_tag(tag, line_num, 'G', i + min_fields, errors)


def validate_o_group(fields: List[str], line_num: int, errors: List[str], defined_ids: Set[str]):
    # O <oid:opt_id> <ref>([ ]<ref>)* <tag>*
    # The spec shows space-separated refs, but GFA is tab-delimited. Assuming refs are in one field, space-separated.
    min_fields = 3
    if len(fields) < min_fields:
        errors.append(f"L{line_num} (O): Expected at least {min_fields} fields, got {len(fields)}.")
        return

    # Field 1: oid (Optional Ordered Group ID)
    oid = fields[1]
    if _validate_opt_id(oid, line_num, 'O', 1, errors, id_type="group ID"):
        if oid != '*' and oid in defined_ids:
            errors.append(f"L{line_num} (O): Duplicate ID definition '{oid}'.")
        elif oid != '*':
            defined_ids.add(oid)

    # Field 2: refs (Space-separated list of references)
    refs_str = fields[2]
    if not refs_str:
         errors.append(f"L{line_num} (O): Reference list (field 2) cannot be empty.")
    else:
        refs = refs_str.split(' ')
        for j, ref_val in enumerate(refs):
            if not _validate_ref(ref_val, line_num, 'O', 2, errors, ref_type=f"reference {j+1}"):
                # Error added by helper
                pass
            # else: # Optional semantic check
            #    base_id = ref_val[:-1]
            #    if base_id not in defined_ids: errors.append(...)

    # Fields 3+: tags
    for i, tag in enumerate(fields[min_fields:]):
        _validate_tag(tag, line_num, 'O', i + min_fields, errors)


def validate_u_group(fields: List[str], line_num: int, errors: List[str], defined_ids: Set[str]):
    # U <uid:opt_id> <id>([ ]<id>)* <tag>*
    # Assuming space-separated IDs in one field.
    min_fields = 3
    if len(fields) < min_fields:
        errors.append(f"L{line_num} (U): Expected at least {min_fields} fields, got {len(fields)}.")
        return

    # Field 1: uid (Optional Unordered Group ID)
    uid = fields[1]
    if _validate_opt_id(uid, line_num, 'U', 1, errors, id_type="group ID"):
        if uid != '*' and uid in defined_ids:
            errors.append(f"L{line_num} (U): Duplicate ID definition '{uid}'.")
        elif uid != '*':
            defined_ids.add(uid)

    # Field 2: ids (Space-separated list of IDs)
    ids_str = fields[2]
    if not ids_str:
         errors.append(f"L{line_num} (U): ID list (field 2) cannot be empty.")
    else:
        ids_in_group = ids_str.split(' ')
        for j, id_val in enumerate(ids_in_group):
             if not _validate_id(id_val, line_num, 'U', 2, errors, id_type=f"ID {j+1}"):
                 # Error added by helper
                 pass
             # else: # Optional semantic check
             #    if id_val not in defined_ids: errors.append(...)

    # Fields 3+: tags
    for i, tag in enumerate(fields[min_fields:]):
        _validate_tag(tag, line_num, 'U', i + min_fields, errors)


# --- Main Validation Function ---

def validate_gfa2(filepath: str) -> List[str]:
    """
    Validates a GFA2 file according to the specification.

    Args:
        filepath: Path to the GFA2 file.

    Returns:
        A list of error messages. Returns an empty list if validation succeeds.
    """
    errors: List[str] = []
    defined_ids: Set[str] = set() # Namespace for S, E, G, O, U IDs
    line_num = 0
    has_fatal_error = False # Flag for encoding/file errors

    try:
        with open(filepath, 'rb') as fb: # Read as bytes first for ASCII check
            byte_content = fb.read()

        # 1. Check for non-ASCII characters
        for i, byte in enumerate(byte_content):
            if byte > 127:
                # Try to decode for context, but might fail if badly encoded
                try:
                    bad_char_context = byte_content[max(0, i-10):min(len(byte_content), i+10)].decode('utf-8', errors='ignore')
                except Exception:
                    bad_char_context = "(cannot decode context)"
                errors.append(f"File Error: Non-ASCII character (byte value {byte}) found near byte offset {i}. Context: '{bad_char_context}'. GFA2 must use only UTF-8 codepoints <= 127.")
                # Decide if this is fatal or continue checking structure
                # has_fatal_error = True # Treat as fatal for strictness
                # break # Stop checking bytes after first error

        # 2. Decode as UTF-8 (strict) - this also ensures it's valid UTF-8
        try:
             file_content = byte_content.decode('utf-8', errors='strict')
        except UnicodeDecodeError as e:
             errors.append(f"File Error: File is not valid UTF-8: {e}")
             has_fatal_error = True


        # 3. Line-by-line structural validation
        if not has_fatal_error:
            lines = file_content.splitlines() # Split into lines
            for line_content in lines:
                line_num += 1
                line = line_content.strip() # Remove leading/trailing whitespace

                if not line: # Skip empty lines
                    continue

                # Check for leading whitespace (already handled by strip, but good practice)
                # if line_content != line_content.lstrip():
                #    errors.append(f"L{line_num}: Line has leading whitespace.")

                record_type = line[0]
                fields = line.split('\t')

                # Check for empty fields (consecutive tabs)
                if any(f == "" for f in fields[1:]): # Check fields after record type
                     # Find first empty field index for better error message
                     first_empty_idx = -1
                     for idx, fld in enumerate(fields):
                          if idx > 0 and fld == "":
                               first_empty_idx = idx
                               break
                     errors.append(f"L{line_num} ({record_type}): Line contains empty field(s) (consecutive tabs) starting at field index {first_empty_idx}.")
                     # Continue validation as best as possible

                if record_type == 'H':
                    validate_header(fields, line_num, errors)
                elif record_type == 'S':
                    validate_segment(fields, line_num, errors, defined_ids)
                elif record_type == 'F':
                    validate_fragment(fields, line_num, errors, defined_ids)
                elif record_type == 'E':
                    validate_edge(fields, line_num, errors, defined_ids)
                elif record_type == 'G':
                    validate_gap(fields, line_num, errors, defined_ids)
                elif record_type == 'O':
                    validate_o_group(fields, line_num, errors, defined_ids)
                elif record_type == 'U':
                    validate_u_group(fields, line_num, errors, defined_ids)
                else:
                    # Ignore lines not starting with recognized codes, as per spec
                    # Optionally add a warning/info message here if desired
                    # print(f"Info L{line_num}: Ignoring line starting with unrecognized code '{record_type}'.")
                    pass

    except FileNotFoundError:
        errors.append(f"File Error: File not found: {filepath}")
        has_fatal_error = True
    except Exception as e:
        errors.append(f"Unexpected Error: An unexpected error occurred during validation: {e}")
        has_fatal_error = True

    # Optional: Add semantic checks here, e.g., all references exist in defined_ids
    # This requires collecting all references during the first pass.

    return errors


# --- Main Execution Block ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate GFA2 file syntax according to the GFA2 specification.",
        epilog="Checks for UTF-8 encoding, ASCII range, record structure, field formats, and duplicate IDs."
    )
    parser.add_argument("gfa2_file", help="Path to the GFA2 file to validate.")
    args = parser.parse_args()

    validation_errors = validate_gfa2(args.gfa2_file)

    if validation_errors:
        print(f"Validation Failed for '{args.gfa2_file}':", file=sys.stderr)
        for error in validation_errors:
            print(f"- {error}", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"Validation Successful: '{args.gfa2_file}' appears to be valid GFA2.", file=sys.stdout)
        sys.exit(0)
