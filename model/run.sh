./${1}.exe ${1}_outs_${2} 1 ${2}
cp specs.json ${1}_outs_${2}/specs.json
../pipeline ${1}_outs_${2}