import sys, aud
cc = aud.cc
if __name__ == "__main__":
    for n in (16, 32, 48, 64, 96):
        try:
            lower = cc.certify("G", n, 8, log=lambda *a: None)
            print(f"G n={n}: certified {lower!r}  (z-plane BIE value quoted in NOTE.md: 1.5237262488052) {'ok' if lower <= 1.5237262488052 else '*** ABOVE ***'}", flush=True)
        except Exception as e:
            print(f"G n={n}: run failed: {e!r}", flush=True)
