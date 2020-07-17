#!/usr/bin/env bash

module-string () {
    raw=$1;
    number=$(printf "%04d" "`echo $raw | cut -d '-' -f 2`")
    echo $raw | sed "s/-.*-/-${number}-/"
}

expected="refine/1.9.2"
raw="refine/1.9.2"
output=$(module-string ${raw})
if [[ "${expected}" == "${output}" ]]; then
    echo pass ${output}
else
    echo "${expected} not ${output}"
    exit 1
fi

expected="refine/1.9.2-0036-g8cb3e9573"
raw="refine/1.9.2-36-g8cb3e9573"
output=$(module-string ${raw})
if [[ "${expected}" == "${output}" ]]; then
    echo pass ${output}
else
    echo "${expected} not ${output}"
    exit 1
fi


