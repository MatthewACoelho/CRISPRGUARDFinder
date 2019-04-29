#!/bin/sh

param_file=./params.nf
submit=no
force=no

fatal()
{
	echo 1>&2 Fatal $*
	exit 1
}

real_path()
{
	local path=`eval echo "$1"`
	local dir=$(dirname "$path")
	if [ -d "$path" ]
	then
		echo $(cd "$path"; pwd)
	elif [ -d "$dir" ]
	then
		echo $(cd "$dir"; pwd)/$(basename "$path")
	else
		echo "$path"
	fi
}

help()
{
	/bin/cat <<END
Usage:
    $0 --help
    $0 [--param_file F] [--submit] [--force]
Options:
    --help            display this
    --force           [optional] re-run everything even if it appears up to date
                      default: no
    --submit          [optional] submit jobs for each sample and finalisation to slurm
                      default: run commands locally
    --param_file F    [optional] file of shell parameters
                      default: $param_file
END
	exit
}

while [ $# -gt 0 ]
do
	case $1 in
	--help)
		help
		;;
	--submit)
		submit=yes
		;;
	--force)
		force=yes
		;;
	-p|--param_file)
		param_file=$( real_path $2 ); shift
		;;
	*)
		fatal "unknown argument: $1"
	esac
	shift
done

[ -d "$guard_root" ] || fatal "No guard_root"
[ -d "$guard_root/bin" ] || fatal "No guard_root/bin"
[ -d "$guard_root/data" ] || fatal "No guard_root/data"
[ -e "$param_file" ] || fatal "Parameter file '$param_file' does not exist"
nf_script=$guard_root/bin/find_guards.nf
[ -e "$nf_script" ] || fatal "Nextflow script file '$nf_script' does not exist"

nf_opts="-C $param_file"
if [ "$submit" = "yes" ]
then
	nf_opts="$nf_opts -C $guard_root/bin/slurm.nf"
else
	nf_opts="$nf_opts -C $guard_root/bin/local.nf"
fi

nf_run_opts=
if [ "$force" = "no" ]
then
	nf_run_opts="$nf_run_opts -resume"
fi

nextflow $nf_opts run $nf_run_opts $nf_script
