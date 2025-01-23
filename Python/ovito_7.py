import os
import gzip
import shutil
from ovito.io import import_file, export_file

def process_files(input_files, output_base_dir='processed_files_4', frame_interval=1000):
    # Ensure the output directory exists
    if not os.path.exists(output_base_dir):
        os.makedirs(output_base_dir)

    # Process each input file
    for infile in input_files:
        # Extract file name without extension (simplified)
        filename = os.path.splitext(os.path.basename(infile))[0]

        # Define output file path
        output_file = os.path.join(output_base_dir, f"{filename}_reduced.dump")

        # Ensure the output file is empty before starting
        with open(output_file, 'w') as f:
            pass

        # Decompress the input file if it's gzipped
        decompressed_file = infile
        if infile.endswith('.gz'):
            decompressed_file = os.path.join(output_base_dir, f"{filename}.lammps")
            with gzip.open(infile, 'rb') as f_in:
                with open(decompressed_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        try:
            # Create a new pipeline for each file
            pipeline = import_file(decompressed_file)
            print(f"Processing {infile}")

            # Export every 1000th frame (combined arguments)
            export_file(pipeline, output_file, format="lammps/dump", multiple_frames=True,
                        columns=["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"],
                        every_nth_frame=frame_interval)
        finally:
            # Delete the decompressed file if it was created
            if infile.endswith('.gz'):
                os.remove(decompressed_file)

# List of input files
#input_files = ['dump_nvt_523_Na10Si2P4O19.out.gz', 'dump_nvt_523_Na12Si6P2O23.out.gz'] 
#input_files = ['dump_nvt_773_Na10Si3P2O14.out.gz', 'dump_nvt_773_Na6Si3P4O19.out.gz']
#input_files = ['dump_nvt_523_Na10Si3P2O14.out.gz', 'dump_nvt_523_Na6Si3P4O19.out.gz']
#input_files = ['dump_nvt_773_Na12Si3P2O17.out.gz', 'dump_nvt_523_Na12Si3P2O17.out.gz']
input_files = ['dump_nvt_773_Na10Si2P4O19.out.gz','dump_nvt_773_Na12Si6P2O23.out.gz']

# Call the function to process files
process_files(input_files, frame_interval=500)


