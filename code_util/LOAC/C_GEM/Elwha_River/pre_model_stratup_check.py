import os

def check_output_writability(file_list):
    """
    Checks if all output files can be opened for writing.
    Does not delete or modify any files.

    Parameters:
        file_list (list of str): List of output filenames.

    Returns:
        bool: True if all files are writable, False otherwise.
    """
    failed = []

    print("Checking writability for output files...\n")

    for fname in file_list:
        full_path = os.path.abspath(fname)

        try:
            with open(full_path, "a"):
                print(f"\033[92mWritable:\033[0m {full_path}")
        except (PermissionError, IOError) as e:
            print(f"ERROR: Cannot write to '{full_path}': {e}")
            failed.append(full_path)

    if failed:
        print("\n\033[31m The following files are not writable. Please close them or check permissions:\033[0m")
        for f in failed:
            print(f"  - {f}")
        return False

    print("\n\033[92mAll output files are writable.\033[0m")
    return True
