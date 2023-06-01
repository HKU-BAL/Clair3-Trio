ori_dir="/autofs/bal36/jhsu/r10/output/4_build_tensors/build/ALL_1356/bins/"

new_dir="/autofs/bal36/jhsu/r10/output/4_build_tensors/build/merged_bins"

mkdir -p ${new_dir}

_dep=10
_chr=(1 2 3 4 5 6 7 8 9)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_chr=(10 11 12 13 14 15)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 2 3 10`

_chr=(16 17 18 19 21 22)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 3 3 10`



_dep=30
_chr=(10 11 12 13 14 15)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_chr=(16 17 18 19 21 22)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 2 3 10`

_chr=(1 2 3 4 5 6 7 8 9)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 3 3 10`



_dep=50
_chr=(10 11 12 13 14 15)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_chr=(1 2 3 4 5 6 7 8 9)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 2 3 10`

_chr=(16 17 18 19 21 22)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 3 3 10`



_dep=65
_chr=(16 17 18 7 8 9)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`


ori_dir="/autofs/bal36/jhsu/r10/output/4_build_tensors/build/SUB_1356/bins/"
new_dir="/autofs/bal36/jhsu/r10/output/4_build_tensors/build/merged_bins"


_dep=SUB_31
_chr=(1 2 3 4 5 6 7 8 9)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_dep=SUB_51
_chr=(10 11 12 13 14 15)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`

_dep=SUB_53
_chr=(16 17 18 19 21 22)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 3 3 10`

_dep=SUB_61
_chr=(1 2 3 4)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 2 3 10`

_dep=SUB_63
_chr=(16 17 18)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 2 3 10`

_dep=SUB_6550
_chr=(10 11 12)
parallel "ln -s ${ori_dir}/HG002_TRIO_${_dep}_{1}_{2} ${new_dir}/HG002_TRIO_${_dep}_{1}_{2}" ::: ${_chr[@]} ::: `seq 1 3 10`
