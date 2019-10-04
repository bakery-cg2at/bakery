for d in p7abx_*[0-9]; do
    oldpwd=$PWD
    cd $d
    start_backmapping.py @params &> ${oldpwd}/log_backmapping_${d}
    cd $oldpwd
done
