import subprocess
import os


script_dir = os.path.dirname(os.path.abspath(__file__))

def test_minimal_input():
    result = subprocess.run(
        ["longbow", "-i", f"{script_dir}/data/sample7_R10D0SUP.fastq"],
        capture_output = True,
        text = True             
    )
    
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"
    
    result_dict = eval(result.stdout)
    assert result_dict["Sample"] == "sample7_R10D0SUP.fastq"
    assert result_dict["Flowcell"] == "R10"
    assert result_dict["Software"] == "dorado"
    assert result_dict["Mode"] == "SUP"


def test_stdout():
    result = subprocess.run(
        ["longbow", "-i", f"{script_dir}/data/sample7_R10D0SUP.fastq", "--stdout"],
        capture_output = True,
        text = True
    )

    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    result_dict = eval(result.stdout)
    assert result_dict["Sample"] == "sample7_R10D0SUP.fastq"
    assert result_dict["Flowcell"] == "R10"
    assert result_dict["Software"] == "dorado"
    assert result_dict["Mode"] == "SUP"


def test_output_to_json():
    result = subprocess.run(
        ["longbow", "-i", f"{script_dir}/data/sample7_R10D0SUP.fastq", "-o", f"{script_dir}/sample7_R10D0SUP.fastq.longbow.json"],
        capture_output = True,
        text = True
    )

    assert result.returncode == 0, f"Command failed with error: {result.stderr}"



def test_output_to_json_and_stdout():
    result = subprocess.run(
        ["longbow", "-i", f"{script_dir}/data/sample7_R10D0SUP.fastq", "-o", f"{script_dir}/sample7_R10D0SUP.fastq.longbow.json", "--stdout"],
        capture_output = True,
        text = True
    )

    assert result.returncode == 0, f"Command failed with error: {result.stderr}"
    result_dict = eval(result.stdout)
    assert result_dict["Sample"] == "sample7_R10D0SUP.fastq"
    assert result_dict["Flowcell"] == "R10"
    assert result_dict["Software"] == "dorado"
    assert result_dict["Mode"] == "SUP"


