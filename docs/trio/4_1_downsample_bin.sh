# this script is used to downsample the generated bin files for model training
# this script is optional, but is recommand for efficient model training
# we downsample the generated bin files by selecting partial bin files into merged_bins folder,
# which will be used for model training.


# original bin dataset
ori_dir="/autofs/bal36/jhsu/r10/output/4_build_tensors/build/ALL_1356/bins/"

# downsampled bin dataset
new_dir="XXXXXX/output/4_build_tensors/build/merged_bins"
mkdir -p ${new_dir}

_dep="10_10_10"
_chr=(1 2 3 4 5 6 7 8 9)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_chr=(10 11 12 13 14 15)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 2 3 10`

_chr=(16 17 18 19 21 22)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 3 3 10`



_dep="30_30_30"
_chr=(10 11 12 13 14 15)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_chr=(16 17 18 19 21 22)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 2 3 10`

_chr=(1 2 3 4 5 6 7 8 9)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 3 3 10`



_dep="60_60_60"
_chr=(10 11 12 13 14 15)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_chr=(1 2 3 4 5 6 7 8 9)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 2 3 10`

_chr=(16 17 18 19 21 22)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 3 3 10`



_dep="80_80_80"
_chr=(16 17 18 7 8 9)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`



_dep="30_10_10"
_chr=(1 2 3 4 5 6 7 8 9)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_dep="60_10_10"
_chr=(10 11 12 13 14 15)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_dep="80_10_10"
_chr=(10 11 12)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_dep="60_30_30"
_chr=(16 17 18 19 21 22)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 3 3 10`

_dep="80_30_30"
_chr=(1 2 3 4)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 2 3 10`

_dep="80_60_60"
_chr=(16 17 18)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 2 3 10`

