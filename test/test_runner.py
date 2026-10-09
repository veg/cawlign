import subprocess
import sys
import os

def compare_fasta(file1, file2):
    def read_fasta(file_path):
        if not os.path.exists(file_path):
            return {}
        with open(file_path, 'r') as f:
            lines = f.readlines()
        sequences = {}
        current_name = None
        current_seq = []
        for line in lines:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            sequences[current_name] = "".join(current_seq)
        return sequences

    try:
        seqs1 = read_fasta(file1)
        seqs2 = read_fasta(file2)
    except Exception as e:
        return False, f"Error reading files: {e}"

    if seqs1.keys() != seqs2.keys():
        return False, f"Different sequence names: {list(seqs1.keys())} vs {list(seqs2.keys())}"

    for name in seqs1:
        if seqs1[name] != seqs2[name]:
            return False, f"Sequence mismatch for {name}"

    return True, "Match"

def main():
    args = sys.argv[1:]
    expect_fail = False
    expect_error = None
    
    while args and args[0].startswith("--"):
        arg = args.pop(0)
        if arg == "--expect-fail":
            expect_fail = True
        elif arg == "--expect-error":
            expect_error = args.pop(0)
        else:
            print(f"Unknown test runner flag: {arg}")
            sys.exit(1)

    if len(args) < 2:
        print("Usage: test_runner.py [--expect-fail] [--expect-error MSG] <baseline_file> <command...>")
        sys.exit(1)

    baseline_file = args[0]
    command = args[1:]

    temp_output = "temp_test_output.fa"
    
    try:
        with open(temp_output, 'w') as f:
            result = subprocess.run(command, stdout=f, stderr=subprocess.PIPE, text=True)
        
        # print(f"DEBUG: Command: {' '.join(command)}")
        # print(f"DEBUG: Return code: {result.returncode}")
        # print(f"DEBUG: Stderr: {result.stderr}")

        if expect_fail:
            if result.returncode == 0:
                print(f"Test failed: Expected command to fail, but it succeeded with code 0.")
                print(f"Command was: {' '.join(command)}")
                sys.exit(1)
            if expect_error and expect_error not in result.stderr:
                print(f"Test failed: Could not find expected error message '{expect_error}' in stderr.")
                print(f"Actual stderr: {result.stderr}")
                sys.exit(1)
            print("Test passed (Expected Failure)!")
            sys.exit(0)

        if result.returncode != 0:
            print(f"Command failed unexpectedly with exit code {result.returncode}")
            print(f"Stderr: {result.stderr}")
            sys.exit(1)

        match, message = compare_fasta(baseline_file, temp_output)
        if not match:
            print(f"Test failed: {message}")
            sys.exit(1)
        else:
            print("Test passed!")
            sys.exit(0)
    finally:
        if os.path.exists(temp_output):
            os.remove(temp_output)

if __name__ == "__main__":
    main()
