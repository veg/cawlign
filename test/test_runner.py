import subprocess
import sys
import os
import tempfile

def compare_fasta(file1, file2, allow_empty=False):
    def read_fasta(file_path):
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        with open(file_path, 'r') as f:
            lines = f.readlines()
        sequences = []
        current_name = None
        current_seq = []
        for line in lines:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if current_name is not None:
                    sequences.append((current_name, "".join(current_seq)))
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name is not None:
            sequences.append((current_name, "".join(current_seq)))
        return sequences

    try:
        seqs1 = read_fasta(file1)
        seqs2 = read_fasta(file2)
    except Exception as e:
        return False, f"Error reading files: {e}"

    if len(seqs1) == 0 and len(seqs2) == 0:
        if allow_empty:
            return True, "Match (both empty)"
        return False, "Both input files are empty!"

    if len(seqs1) != len(seqs2):
        return False, f"Different record counts: {len(seqs1)} vs {len(seqs2)}"

    from collections import defaultdict
    dict1 = defaultdict(list)
    dict2 = defaultdict(list)
    for name, seq in seqs1:
        dict1[name].append(seq)
    for name, seq in seqs2:
        dict2[name].append(seq)

    if set(dict1.keys()) != set(dict2.keys()):
        diff1 = set(dict1.keys()) - set(dict2.keys())
        diff2 = set(dict2.keys()) - set(dict1.keys())
        return False, f"Record name mismatch. Only in baseline: {diff1}, only in output: {diff2}"

    for name in dict1:
        s1_list = sorted(dict1[name])
        s2_list = sorted(dict2[name])
        if s1_list != s2_list:
            return False, f"Sequence mismatch for record '{name}'"

    return True, "Match"

def main():
    args = sys.argv[1:]
    expect_fail = False
    expect_error = None
    check_empty_stderr = False
    allow_empty = False
    
    while args and args[0].startswith("--"):
        arg = args.pop(0)
        if arg == "--expect-fail":
            expect_fail = True
        elif arg == "--expect-error":
            expect_error = args.pop(0)
        elif arg == "--check-empty-stderr":
            check_empty_stderr = True
        elif arg == "--allow-empty":
            allow_empty = True
        else:
            print(f"Unknown test runner flag: {arg}")
            sys.exit(1)

    if len(args) < 2:
        print("Usage: test_runner.py [--expect-fail] [--expect-error MSG] [--check-empty-stderr] [--allow-empty] <baseline_file> <command...>")
        sys.exit(1)

    baseline_file = args[0]
    command = args[1:]

    with tempfile.NamedTemporaryFile(suffix='.fa', delete=False) as tmp:
        temp_output = tmp.name
    
    try:
        with open(temp_output, 'w') as f:
            result = subprocess.run(command, stdout=f, stderr=subprocess.PIPE, text=True)

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

        if check_empty_stderr and result.stderr.strip():
            print(f"Test failed: Expected empty stderr in quiet mode, got: {result.stderr}")
            sys.exit(1)

        match, message = compare_fasta(baseline_file, temp_output, allow_empty=allow_empty)
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
