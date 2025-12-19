import filecmp
import subprocess


def test_parse_timestamps_stdout(path_timestamp_a1, tmp_path):
    t1 = tmp_path / "intermediate1.ts"
    t2 = tmp_path / "intermediate2.ts"
    t3 = tmp_path / "intermediate3.ts"
    t4 = tmp_path / "intermediate4.ts"

    subprocess.check_call(f"parse-timestamps -X {path_timestamp_a1} - >{t1}", shell=True)
    subprocess.check_call(f"parse-timestamps -X {path_timestamp_a1} {t2}", shell=True)
    subprocess.check_call(f"cat {path_timestamp_a1} | parse-timestamps -X - {t3}", shell=True)
    subprocess.check_call(f"cat {path_timestamp_a1} | parse-timestamps -X - - >{t4}", shell=True)

    assert filecmp.cmp(t1, t2)
    assert filecmp.cmp(t2, t3)
    assert filecmp.cmp(t3, t4)

def test_parse_timestamps_stdout_a0(path_timestamp_a0, tmp_path):
    t5 = tmp_path / "intermediate5.ts"
    t6 = tmp_path / "intermediate6.ts"
    subprocess.check_call(f"parse-timestamps -X -A0 -a0 {path_timestamp_a0} - >{t5}", shell=True)
    subprocess.check_call(f"parse-timestamps -X -A0 -a0 {path_timestamp_a0} {t6}", shell=True)
    assert filecmp.cmp(t5, t6)

def test_parse_timestamps_stdout_a2(path_timestamp_a2, tmp_path):
    t7 = tmp_path / "intermediate7.ts"
    t8 = tmp_path / "intermediate8.ts"
    subprocess.check_call(f"parse-timestamps -X -A2 -a2 {path_timestamp_a2} - >{t7}", shell=True)
    subprocess.check_call(f"parse-timestamps -X -A2 -a2 {path_timestamp_a2} {t8}", shell=True)
    assert filecmp.cmp(t7, t8)
