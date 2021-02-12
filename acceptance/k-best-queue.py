#!/usr/bin/env python

print('K3-debug')
exit()

from subprocess import check_output
import re

if __name__ == "__main__":
    qstat_output = check_output("qstat -a | grep Free | egrep -v 'hmem|KNL|K2'", shell=True)
    free_nodes = {}
    found_queue = False
    for line in qstat_output.split("\n"):
        m = re.match("^\s*(\d+).*--\s+(K\d)\s+\(ALL.*Free", line)
        if m:
            found_queue = True
            free_nodes[m.groups(0)[1]] = {"nodes available": int(m.groups(0)[0])}

    if not found_queue:
        print "K2a-fun3d"
        exit(0)

    for q in free_nodes:
        queue_output = check_output("qstat -q | grep %s-standard | grep -v 512" % q, shell=True)
        free_nodes[q]["selection weight"] = free_nodes[q]["nodes available"] - int(queue_output.split()[6])

        # Try to get K4 more often
        if q == "K4":
            if free_nodes[q]["nodes available"] > 10:
                print "K4-debug"
                exit(0)

    cores_count = {"K4": 40, "K3": 16, "K2": 12, "K2a": 12}

    highest_weight = -99999
    best_queue_to_submit = ""
    for q in free_nodes:
        w = free_nodes[q]["selection weight"]
        w = max(w, w * cores_count[q])
        if w > highest_weight:
            highest_weight = w
            best_queue_to_submit = q + "-debug"

    if highest_weight <= -20:
        best_queue_to_submit = "K2a-fun3d"

    print best_queue_to_submit
