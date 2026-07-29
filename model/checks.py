class Checks:
    def __init__(self):
        self.fails = 0

    def __call__(self, name, ok, detail=""):
        self.fails += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))

    def report(self):
        print(f"\n{'ALL PASS' if not self.fails else f'{self.fails} FAILED'}")
        raise SystemExit(1 if self.fails else 0)
