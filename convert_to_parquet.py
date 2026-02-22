# convert_to_parquet.py
import pyarrow as pa
import pyarrow.csv as pv
import pyarrow.parquet as pq
import argparse
import time
import gzip
from tqdm import tqdm

def get_headers(filepath):
    """Reads just the first line to get the exact column names."""
    print("[Info] Extracting headers to build schema...")
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        header_line = f.readline().strip('\n')
        return header_line.split('\t')

def main():
    parser = argparse.ArgumentParser(description="Convert TSV.gz to Parquet (Low Memory Streaming)")
    parser.add_argument("-i", "--input", required=True, help="Input TSV.gz file")
    parser.add_argument("-o", "--output", required=True, help="Output Parquet file")
    parser.add_argument("-t", "--total", type=int, default=None, 
                        help="Total number of rows in the input file (to calculate precise ETA)")
    args = parser.parse_args()

    print(f"[Start] Converting {args.input} to {args.output} (Streaming Mode)...")
    start_time = time.time()

    # 1. Force all columns to be read as strings to prevent Null inference errors on sparse data
    headers = get_headers(args.input)
    column_types = {col: pa.string() for col in headers}
    
    # 2. Configure options
    parse_options = pv.ParseOptions(delimiter='\t')
    read_options = pv.ReadOptions(block_size=50 * 1024 * 1024) 
    convert_options = pv.ConvertOptions(column_types=column_types)

    print("[Info] Initializing stream. Memory usage will remain flat.")
    
    # 3. Stream and write
    with pv.open_csv(args.input, parse_options=parse_options, read_options=read_options, convert_options=convert_options) as reader:
        writer = None
        
        with tqdm(total=args.total, desc="Converting", unit=" rows") as pbar:
            for batch in reader:
                if writer is None:
                    writer = pq.ParquetWriter(
                        args.output, 
                        batch.schema,
                        compression='snappy',
                        use_dictionary=True  # Compresses repeated strings highly efficiently
                    )
                
                writer.write_batch(batch)
                pbar.update(batch.num_rows)

        if writer:
            writer.close()

    print(f"\n[Done] Conversion finished in {(time.time() - start_time) / 60:.2f} minutes.")

if __name__ == "__main__":
    main()