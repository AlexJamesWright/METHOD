TOL=10e-15

time_format_folder="Final"
vars_folders=["Auxiliary", "Conserved", "Primitive"]
vars_files=["aux", "cons", "prims"]
extension=".dat"

for index in range(len(vars_folders)):
        serial_filename = "/".join(["Serial", "Data", time_format_folder, vars_folders[index], vars_files[index]])
        parallel_filename = "/".join(["ParallelMPI", "Data", time_format_folder, vars_folders[index], vars_files[index]])
        serial_filename = serial_filename+extension
        parallel_filename = parallel_filename+extension
        print("Processing: " + serial_filename + ", " + parallel_filename)

        try:
                with open(serial_filename, 'r') as serial_dat_file:
                        with open(parallel_filename, 'r') as parallel_dat_file:
                                skip_header = 1
                                for serial_line, parallel_line in zip(serial_dat_file, parallel_dat_file):
                                    if skip_header:
                                        skip_header = 0
                                        continue
                                    serial_val = float(serial_line)
                                    parallel_val = float(parallel_line)
                                    if (abs(serial_val-parallel_val) > TOL):
                                        print("\n\n!! Error in {} (val={}), {}, (val={})\n\n".format(serial_filename, serial_val, parallel_filename, parallel_val))  
                                        break

        except IOError:
                print("Could not read file:", filename)


