# MultiscaleModeling

Internal collaboration on multiscale transmission modeling projects

**Models that live outside this repo.**
[Intrahost polio shedding, susceptibility, and immunity model](https://github.com/famulare/cessationStability)

### Profiling

Install Snakeviz to visualize the output file of the profiler.
```bash
pip install snakeviz
```

Then run the profiler.
```bash
python.exe -B -m cProfile -o output_profiler.prof outbreak_local.py
```

Snakeviz will run a server and display the result in the default browser.
```bash
snakeviz output_profiler.prof
```

Another intresting tool for profiling is vprof (https://github.com/nvdv/vprof), it can display the time spent per line.

```bash
pip install vprof
```

Running vprof with parameter h will display a code heatmap in the browser.

```bash
vprof -c h outbreak_local.py
```

vprof adds quite some overhead so it is a good idea to keep the program short.