using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
# where the standard output will be saved
output_dir = projectdir()
o = output_dir*"/IO/o"  #where the standard output of the program is saved
e = output_dir*"/IO/e"  #where the error messages are saved
#(Notice that my code has its file output configured from within)

isdir(e) || mkpath(e)
isdir(o) || mkpath(o)

# Location of the scripts to-be-submitted
scriptdir = output_dir*"scripts/"
scripts = ["mpsb/mpsb_mct.jl", "mpsb/mpsb_lyapunov.jl"]

# Queues:
queue = "rostam.q"

for script in scripts
  name = "$(basename(script))"

  julia_submit = "julia $(joinpath(scriptdir, script)) "
  command = `qsub -N $name -o $o -e $e -q $queue -cwd -b y $julia_submit`

  run(command)
end

# TESTING
# julia_submit = "/usr/nld/julia-0.6.0/bin/julia $(script) 0"
# command = `qsub -o $o -e $e -q $queue -cwd -b y $julia_submit`
# # command = `$julia_submit`
# run(command)
