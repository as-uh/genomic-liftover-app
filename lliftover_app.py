import streamlit as st
import pandas as pd
import numpy as np
import os
import tempfile
import subprocess # For more robust command execution than os.system

# --- Original Helper Functions (slightly adapted for clarity if needed) ---
def create_marker_id(row, chrom_col, pos_col, ref_allele_col, alt_allele_col):
    """Creates a unique marker ID from chromosome, position, and alleles."""
    try:
        # Ensure alleles are strings before sorting, in case they are numeric
        ref_allele = str(row[ref_allele_col])
        alt_allele = str(row[alt_allele_col])
        alleles = sorted([ref_allele, alt_allele])
        return f"chr{row[chrom_col]}:{row[pos_col]}:{alleles[0]}:{alleles[1]}"
    except KeyError as e:
        st.error(f"Column not found: {e}. Please check your column name inputs.")
        st.stop()
    except Exception as e:
        st.error(f"Error creating marker ID for row: {row}\nError: {e}")
        # Optionally, return a placeholder or skip the row
        return f"ERROR_CREATING_ID_chr{row.get(chrom_col, 'NA')}:{row.get(pos_col, 'NA')}"


def convert_chromosome(chrom):
    """Converts chromosome values to integers. Handles 'chr' prefix."""
    if isinstance(chrom, str):
        chrom_val = chrom.replace('chr', '').strip()
        if chrom_val == 'X':
            return 23
        elif chrom_val == 'Y':
            return 24
        elif chrom_val == 'MT' or chrom_val == 'M':
            return 25 # Common convention for mitochondrial
    else: # if it's already numeric
        chrom_val = chrom

    try:
        return int(chrom_val)
    except ValueError:
        return np.nan


def create_bed_file(data, chrom_col, pos_col, marker_id_col, output_bed_path):
    """Creates a BED file for CrossMap."""
    bed = data.copy()
    # Ensure 'marker_id' column exists from create_marker_id
    if marker_id_col not in bed.columns:
        st.error(f"'{marker_id_col}' column not found. Marker ID creation might have failed.")
        st.stop()

    # Convert numeric chromosome back to string for BED 'chrX' format
    # The convert_chromosome function is for converting TO numeric for internal use
    # For BED, we need 'chr1', 'chrX' etc.
    def format_chrom_for_bed(c):
        if c == 23: return 'chrX'
        if c == 24: return 'chrY'
        if c == 25: return 'chrM' # or chrMT depending on chain file
        return f'chr{c}'

    try:
        # Assuming chrom_col contains numeric or X/Y which was converted to numeric
        # If original chrom_col had 'chr' prefix, it needs to be handled before numeric conversion
        # For simplicity, let's assume chrom_col is already prepared as numbers (1-22, 23 for X, 24 for Y)
        # If not, an earlier step should convert it.
        # Let's assume the 'chrom_col' in the input DataFrame is what we need to format.
        bed['bed_chrom'] = bed[chrom_col].apply(lambda x: format_chrom_for_bed(convert_chromosome(x)))

        bed['start'] = bed[pos_col].astype(int) -1 # BED is 0-indexed start
        bed['end'] = bed[pos_col].astype(int)
        bed_df_to_save = bed[['bed_chrom', 'start', 'end', marker_id_col]]

        # Remove rows where bed_chrom could not be determined (e.g. 'chrNA')
        bed_df_to_save = bed_df_to_save[~bed_df_to_save['bed_chrom'].str.contains("NA", na=False)]
        bed_df_to_save = bed_df_to_save.dropna(subset=['bed_chrom', 'start', 'end', marker_id_col])


        if bed_df_to_save.empty:
            st.error("No valid data to write to BED file after formatting. Check chromosome/position columns.")
            st.stop()

        bed_df_to_save.to_csv(output_bed_path, sep='\t', index=False, header=False)
        st.info(f"BED file created at {output_bed_path} with {len(bed_df_to_save)} entries.")

    except KeyError as e:
        st.error(f"Column not found while creating BED file: {e}. Ensure '{chrom_col}', '{pos_col}', and '{marker_id_col}' are correct.")
        st.stop()
    except Exception as e:
        st.error(f"An error occurred during BED file creation: {e}")
        st.stop()


def lift_over_coordinates(data_path, chain_path, output_prefix,
                          chrom_col, pos_col, ref_col, alt_col,
                          sep='\t', compression='infer'):
    """Performs liftover from one genome build to another using CrossMap."""
    try:
        data = pd.read_csv(data_path, sep=sep, compression=compression, low_memory=False)
        st.write(f"Original data loaded: {len(data)} rows.")
        if data.empty:
            st.error("Input data file is empty or could not be read correctly.")
            return None
    except Exception as e:
        st.error(f"Error reading data file: {e}")
        return None

    # Convert chromosome column to a standardized numeric representation early
    # This helps if input has 'chrX' or just 'X' or '23'
    data['numeric_chrom'] = data[chrom_col].apply(convert_chromosome)
    data = data.dropna(subset=['numeric_chrom']) # Remove rows where chromosome couldn't be parsed
    data['numeric_chrom'] = data['numeric_chrom'].astype(int)
    st.write(f"Data after chromosome conversion & filtering: {len(data)} rows.")
    if data.empty:
        st.error("No valid chromosome entries found after initial processing.")
        return None

    # Create marker_id using the standardized numeric chromosome
    # but the marker_id function itself will prepend 'chr' based on the original chrom_col value
    # Let's ensure create_marker_id uses a consistent chromosome representation for the ID
    # So, we pass the original 'chrom_col' for the ID string, but use 'numeric_chrom' for BED
    temp_marker_id_col = '_temp_marker_id'
    data[temp_marker_id_col] = data.apply(
        lambda row: create_marker_id(row, chrom_col, pos_col, ref_col, alt_col),
        axis=1
    )
    st.write("Marker IDs created.")

    # Define BED file paths
    input_bed_path = f"{output_prefix}.bed"
    lifted_bed_path = f"{output_prefix}_lifted.bed" # CrossMap output
    unmapped_bed_path = f"{output_prefix}_unmapped.bed" # CrossMap unmapped output

    # Create BED file using 'numeric_chrom' for consistency in BED chr format
    create_bed_file(data, 'numeric_chrom', pos_col, temp_marker_id_col, input_bed_path)

    # Run CrossMap
    # Ensure CrossMap is in PATH or provide full path
    crossmap_command = [
        "CrossMap", "bed",
        chain_path,
        input_bed_path,
        lifted_bed_path
    ]
    # Adding unmapped output
    # CrossMap.py bed [OPTIONS] CHAIN_FILE INPUT_FILE OUTPUT_FILE
    # If you want unmapped: CrossMap.py bed [OPTIONS] CHAIN_FILE INPUT_FILE OUTPUT_FILE_PREFIX
    # In this case, INPUT_FILE is bed, OUTPUT_FILE_PREFIX will generate .bed and .unmap.bed
    # So the command should be:
    # CrossMap bed chain_file input_bed_path output_prefix (this will create output_prefix.bed and output_prefix.unmap.bed)
    # Let's adjust:
    lifted_bed_path = f"{output_prefix}_lifted.bed" # Explicitly name lifted output
    # CrossMap generates output_prefix.bed, so we rename it to be clear
    # Or, if CrossMap directly takes the output file name:
    # CrossMap bed {chain_file} {input_bed} {output_lifted_bed}
    # The original script had: CrossMap bed {chain_path} {output_path}.bed {output_path}_hg38.bed
    # This means it writes mapped to {output_path}_hg38.bed and unmapped to {output_path}_hg38.bed.unmap
    
    # Let's stick to the original script's naming for CrossMap output for mapped reads
    actual_crossmap_output_bed = f"{output_prefix}_crossmap_output.bed" # This is where CrossMap writes mapped
    unmapped_crossmap_file = f"{actual_crossmap_output_bed}.unmap" # CrossMap appends .unmap

    st.info(f"Running CrossMap: {' '.join(crossmap_command)}")
    st.info(f"Input BED for CrossMap: {input_bed_path}")
    st.info(f"Chain file for CrossMap: {chain_path}")
    st.info(f"Output BED from CrossMap will be: {actual_crossmap_output_bed}")


    # Command: CrossMap bed <chain_file> <input_bed> <output_bed>
    # The third argument is the output file for *mapped* regions.
    # Unmapped regions are typically written to <output_bed>.unmap
    crossmap_command_final = ["CrossMap", "bed", chain_path, input_bed_path, actual_crossmap_output_bed]

    try:
        # Using subprocess.run for better control
        result = subprocess.run(crossmap_command_final, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            st.error(f"CrossMap failed with return code {result.returncode}.")
            st.error(f"CrossMap STDERR:\n{result.stderr}")
            st.error(f"CrossMap STDOUT:\n{result.stdout}")
            # Check if CrossMap exists
            try:
                subprocess.run(["CrossMap", "-h"], capture_output=True, check=True)
            except FileNotFoundError:
                st.error("CrossMap command not found. Please ensure it is installed and in your system's PATH.")
            except subprocess.CalledProcessError:
                st.warning("CrossMap found, but '-h' failed. Check CrossMap installation.")
            return None
        st.success("CrossMap finished.")
        st.text(f"CrossMap STDOUT:\n{result.stdout}")
        if result.stderr:
             st.text(f"CrossMap STDERR:\n{result.stderr}")


    except FileNotFoundError:
        st.error("CrossMap command not found. Please ensure it is installed and in your system's PATH.")
        return None
    except Exception as e:
        st.error(f"An error occurred while running CrossMap: {e}")
        return None

    if not os.path.exists(actual_crossmap_output_bed) or os.path.getsize(actual_crossmap_output_bed) == 0:
        st.warning(f"Lifted BED file '{actual_crossmap_output_bed}' was not created or is empty. "
                   "This might mean no coordinates were successfully lifted over. "
                   "Check CrossMap logs if available, and the unmapped file.")
        # Create an empty DataFrame for lifted if file doesn't exist to avoid crashing later
        lifted = pd.DataFrame(columns=['chrom', 'start', 'end', temp_marker_id_col, 'Pos_b_lifted'])

    else:
        try:
            lifted = pd.read_csv(
                actual_crossmap_output_bed, sep='\t', header=None,
                names=['chrom_lifted_raw', 'start_lifted', 'end_lifted', temp_marker_id_col]
            )
            st.write(f"Lifted data read from BED: {len(lifted)} rows.")
            # Position in BED 'end' column is the 1-based coordinate for single points
            lifted['Pos_b_lifted'] = lifted['end_lifted'].astype('int64')

        except pd.errors.EmptyDataError:
            st.warning(f"Lifted BED file '{actual_crossmap_output_bed}' is empty. No coordinates mapped.")
            lifted = pd.DataFrame(columns=['chrom_lifted_raw', 'start_lifted', 'end_lifted', temp_marker_id_col, 'Pos_b_lifted'])
        except Exception as e:
            st.error(f"Error reading lifted BED file '{actual_crossmap_output_bed}': {e}")
            return None


    if not lifted.empty:
        lifted['chrom_lifted_numeric'] = lifted['chrom_lifted_raw'].apply(convert_chromosome)
        lifted = lifted[lifted['chrom_lifted_numeric'].notnull()]
        if not lifted.empty:
            lifted['chrom_lifted_numeric'] = lifted['chrom_lifted_numeric'].astype('int64')
    else:
        # Add empty columns if lifted is empty to prevent merge errors
        lifted['chrom_lifted_numeric'] = pd.Series(dtype='int64')
        lifted['Pos_b_lifted'] = pd.Series(dtype='int64')
        lifted['chrom_lifted_raw'] = pd.Series(dtype='object')


    st.write(f"Processed lifted data: {len(lifted)} rows.")

    # Merge back with original data
    merged = data.merge(
        lifted[['chrom_lifted_numeric', 'Pos_b_lifted', 'chrom_lifted_raw', temp_marker_id_col]],
        on=temp_marker_id_col,
        how='left'
    )
    # Rename columns for clarity
    merged = merged.rename(columns={
        'chrom_lifted_numeric': 'chrom_lifted',
        'Pos_b_lifted': 'pos_lifted',
        'chrom_lifted_raw': 'chrom_lifted_str' # Keep the chrX style string if needed
    })
    # Remove the temporary numeric chromosome column if it's a duplicate or not needed
    # merged = merged.drop(columns=['numeric_chrom']) # numeric_chrom was from original data
    merged = merged.drop(columns=[temp_marker_id_col]) # Drop the temporary marker ID

    st.write(f"Merged data: {len(merged)} rows.")

    output_csv_path = f"{output_prefix}_lifted_results.csv"
    try:
        merged.to_csv(output_csv_path, sep=sep if sep != '\t' else ',', index=False, compression=compression if compression != "infer" else None)
        st.success(f"Liftover complete. Results saved to {output_csv_path} (for download).")
        return output_csv_path
    except Exception as e:
        st.error(f"Error saving final CSV: {e}")
        return None

# --- Streamlit App UI ---
st.set_page_config(layout="wide")
st.title("Genomic Coordinate Liftover Tool (using CrossMap)")

st.markdown("""
This tool performs a liftover of genomic coordinates from one genome build to another
using a chain file and the CrossMap utility.

**Prerequisite:** `CrossMap` (and its dependency `pysam`) must be installed in the environment
where this Streamlit app is running and be accessible from the command line.
You can typically install it via pip: `pip install CrossMap pysam`
""")

# File Uploads
st.sidebar.header("1. Upload Files")
uploaded_data_file = st.sidebar.file_uploader("Upload your input data file (e.g., .tsv.gz, .csv)", type=['csv', 'tsv', 'gz', 'txt'])
uploaded_chain_file = st.sidebar.file_uploader("Upload your chain file (e.g., .chain.gz)", type=['chain', 'gz'])

# Configuration
st.sidebar.header("2. Configure Parameters")
output_prefix = st.sidebar.text_input("Output file prefix", "liftover_results")

# Dynamically get column names from uploaded data file if available
col_names = []
temp_data_for_cols = None
if uploaded_data_file:
    try:
        # Read only a few lines to get headers
        if uploaded_data_file.name.endswith('.gz'):
            temp_data_for_cols = pd.read_csv(uploaded_data_file, sep=None, engine='python', compression='gzip', nrows=5)
        else:
            temp_data_for_cols = pd.read_csv(uploaded_data_file, sep=None, engine='python', nrows=5)
        uploaded_data_file.seek(0) # Reset file pointer
        if temp_data_for_cols is not None:
            col_names = temp_data_for_cols.columns.tolist()
    except Exception as e:
        st.sidebar.warning(f"Could not pre-read column names: {e}")


chrom_col = st.sidebar.selectbox("Chromosome column name", col_names, index=col_names.index('chrom') if 'chrom' in col_names else (col_names.index('CHR') if 'CHR' in col_names else 0) if col_names else 0)
pos_col = st.sidebar.selectbox("Base pair position column name", col_names, index=col_names.index('pos') if 'pos' in col_names else (col_names.index('POS') if 'POS' in col_names else (col_names.index('BP') if 'BP' in col_names else 0)) if col_names else 0)
ref_col = st.sidebar.selectbox("Reference allele column name", col_names, index=col_names.index('ref') if 'ref' in col_names else (col_names.index('REF') if 'REF' in col_names else 0) if col_names else 0)
alt_col = st.sidebar.selectbox("Alternate allele column name", col_names, index=col_names.index('alt') if 'alt' in col_names else (col_names.index('ALT') if 'ALT' in col_names else 0) if col_names else 0)

input_sep_options = {',': 'Comma (,)', '\t': 'Tab (\\t)', ' ': 'Space ( )'}
input_sep_display = st.sidebar.selectbox("Input file separator", options=list(input_sep_options.keys()), format_func=lambda x: input_sep_options[x])

# Compression for output (input compression is inferred by pandas)
# output_compression = st.sidebar.selectbox("Output file compression", [None, 'gzip'], format_func=lambda x: x if x else "None")
# For simplicity, let's make output CSV uncompressed, or match input if it was gzip
# The lift_over_coordinates function already handles compression='infer' for input and 'None' or 'gzip' for output.
# For the web app, providing an uncompressed CSV for download is often easiest.


if st.sidebar.button("🚀 Run Liftover"):
    if uploaded_data_file and uploaded_chain_file and chrom_col and pos_col and ref_col and alt_col and output_prefix:
        with st.spinner("Processing... This may take a while for large files."):
            # Create temporary directory to store uploaded files and intermediate files
            with tempfile.TemporaryDirectory() as tmpdir:
                data_file_path = os.path.join(tmpdir, uploaded_data_file.name)
                with open(data_file_path, "wb") as f:
                    f.write(uploaded_data_file.getbuffer())

                chain_file_path = os.path.join(tmpdir, uploaded_chain_file.name)
                with open(chain_file_path, "wb") as f:
                    f.write(uploaded_chain_file.getbuffer())

                temp_output_prefix = os.path.join(tmpdir, output_prefix)

                # Determine compression for output based on input data file name
                output_compression_type = 'gzip' if uploaded_data_file.name.endswith('.gz') else None

                final_output_csv_path = lift_over_coordinates(
                    data_path=data_file_path,
                    chain_path=chain_file_path,
                    output_prefix=temp_output_prefix,
                    chrom_col=chrom_col,
                    pos_col=pos_col,
                    ref_col=ref_col,
                    alt_col=alt_col,
                    sep=input_sep_display,
                    compression=output_compression_type # Pass this for saving output
                )

                if final_output_csv_path and os.path.exists(final_output_csv_path):
                    st.success("Liftover process completed!")
                    with open(final_output_csv_path, "rb") as fp:
                        # Determine the final download filename
                        download_filename = f"{output_prefix}_lifted_results.csv"
                        if output_compression_type == 'gzip':
                            download_filename += ".gz"

                        st.download_button(
                            label="Download Lifted Data",
                            data=fp,
                            file_name=download_filename,
                            mime="text/csv" if output_compression_type is None else "application/gzip"
                        )
                    
                    # Optionally, display a preview of the results
                    try:
                        preview_df = pd.read_csv(final_output_csv_path, sep=input_sep_display if input_sep_display != '\t' else ',', 
                                                 compression=output_compression_type)
                        st.dataframe(preview_df.head())
                        
                        # Provide info on unmapped file if it exists
                        unmapped_file_expected_path = f"{temp_output_prefix}_crossmap_output.bed.unmap"
                        if os.path.exists(unmapped_file_expected_path) and os.path.getsize(unmapped_file_expected_path) > 0:
                            st.info(f"An unmapped regions file was also generated by CrossMap: "
                                    f"'{os.path.basename(unmapped_file_expected_path)}'. "
                                    f"This is not included in the download but indicates regions that couldn't be lifted.")
                        elif os.path.exists(unmapped_file_expected_path):
                             st.info("Unmapped regions file was generated by CrossMap but is empty.")


                    except Exception as e:
                        st.warning(f"Could not display preview of results: {e}")

                elif final_output_csv_path: # Path returned but file not found
                     st.error(f"Output file {final_output_csv_path} was expected but not found. Please check logs.")
                else: # None returned, error already shown by function
                    st.error("Liftover process failed. Please check the error messages above.")
    else:
        st.sidebar.error("Please fill in all fields and upload necessary files.")

st.markdown("---")
st.markdown("by Alok Singh - 2025")

st.markdown("If you encounter any issues, please report them on the GitHub repository. https://github.com/as-uh/genomic-liftover-app.git")
st.markdown("**Disclaimer:** This tool is provided as-is. Ensure you have the necessary permissions to use the data and tools.")
st.markdown("**Note:** This app is designed for educational and research purposes. "
             "Always validate results with appropriate bioinformatics pipelines.")
```
# End of the Streamlit app code
```