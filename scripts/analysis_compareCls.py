ranks=["eukaryota", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
diff = {"ssu_lsu": {}, "ssu_its": {}, "its_lsu":{}}
size = {}
allCls = {}
com = {}
stat = {}
fungi = {}
for line in open(snakemake.input.cls):
    oId, tSize, ssuCls, itsCls, lsuCls = line.strip().split("\t")
    size[oId] = int(tSize)
    cls = {}
    maxLen=0
    if ssuCls == "unknown":
        cls["ssu"] = None
    else:
        cls["ssu"] = ssuCls.split(";")
        if len(cls["ssu"]) > maxLen:
            maxLen=len(cls["ssu"])
    if lsuCls == "unknown":
        cls["lsu"] = None
    else:
        cls["lsu"] = [lsuEntry.replace(" ", "_") for lsuEntry in lsuCls.split(";")]
        if len(cls["lsu"]) > maxLen:
            maxLen=len(cls["lsu"])
    if itsCls == "unknown":
        cls["its"] = None
    else:
        cls["its"] = ["Eukaryota"] + [c.split("__", 1)[-1] for c in itsCls.split(";")]
        if len(cls["its"]) > maxLen:
            maxLen=len(cls["its"])
    fungi[oId] = any([cls["ssu"] and len(cls["ssu"]) > 1 and cls["ssu"][1] == "Fungi", 
                      cls["its"] and len(cls["its"]) > 1 and cls["its"][1] == "Fungi", 
                      cls["lsu"] and len(cls["lsu"]) > 1 and cls["lsu"][1] == "Fungi"])
    allCls[oId] = cls
    com[oId] = []
    stat[oId] = {}
    if not cls["ssu"] is None:
        stat[oId]["ssu"] = {"its": ["NA"]*8,
                            "lsu": ["NA"]*8}
    if not cls["its"] is None:
        stat[oId]["its"] = {"ssu": ["NA"]*8,
                            "lsu": ["NA"]*8}
    if not cls["lsu"] is None:
        stat[oId]["lsu"] = {"ssu": ["NA"]*8,
                            "its": ["NA"]*8}
    for lvl in range(maxLen):
        tCls = [None, None, None]
        if not cls["ssu"] is None and len(cls["ssu"]) > lvl:
            tCls[0] = cls["ssu"][lvl]
        if not cls["its"] is None and len(cls["its"]) > lvl:
            tCls[1] = cls["its"][lvl]
        if not cls["lsu"] is None and len(cls["lsu"]) > lvl:
            tCls[2] = cls["lsu"][lvl]
        if not tCls[0] is None and (tCls[0] == tCls[1] or tCls[0] == tCls[2]):
            com[oId].append(tCls[0])
        elif not tCls[1] is None and tCls[1] == tCls[2]:
            com[oId].append(tCls[1])
        else:
            com[oId].append(None)
        for i, mrk in [(0, "ssu"), (1, "its"), (2, "lsu")]:
            if not tCls[i] is None:
                (o1, oMrk1), (o2, oMrk2) = set([(0, "ssu"), (1, "its"), (2, "lsu")]) - set([(i,mrk)])
                if tCls[o1] is None:
                    stat[oId][mrk][oMrk1][lvl] = "blank"
                else:
                    if tCls[i] == tCls[o1]:
                        stat[oId][mrk][oMrk1][lvl] = "same"
                    else:
                        stat[oId][mrk][oMrk1][lvl] = "diff"
                if tCls[o2] is None:
                    stat[oId][mrk][oMrk2][lvl] = "blank"
                else:
                    if tCls[i] == tCls[o2]:
                        stat[oId][mrk][oMrk2][lvl] = "same"
                    else:
                        stat[oId][mrk][oMrk2][lvl] = "diff"
            
    for mrk1, mrk2 in [("ssu","lsu"), ("ssu", "its"), ("its", "lsu")]:
        if cls[mrk1] is None or cls[mrk2] is None:
            continue
        for r, rank in enumerate(ranks):
            if r >= len(cls[mrk1]) or r >= len(cls[mrk2]):
                break
            if cls[mrk1][r] != cls[mrk2][r]:
                if r == 0:
                    parent="root"
                else:
                    parent=cls[mrk1][r-1]
                if r>2:
                    pphylum=cls[mrk1][2]
                else:
                    pphylum="NA"
                if rank not in diff["%s_%s" % (mrk1, mrk2)]:
                    diff["%s_%s" % (mrk1, mrk2)][rank] = []
                diff["%s_%s" % (mrk1, mrk2)][rank].append((oId, cls[mrk1][r], cls[mrk2][r], parent, pphylum))
                break

with open(snakemake.output.diff, "w") as out:
    for comp in diff.keys():
        mrk1, mrk2 = comp.split("_")
        for rank, diffData in diff[comp].items():
            for oId, m1Cls, m2Cls, parent, phylum in diffData:
                out.write("%s\t%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s\t%s\n" % (comp, mrk1, mrk2, rank, oId, size[oId], m1Cls, m2Cls, parent, phylum))
with open(snakemake.output.comp, "w") as out:
    for oId in com.keys():
        out.write("%s\t%s;" % (oId, ";".join([str(c) for c in com[oId]])))
        for marker in ["ssu", "its", "lsu"]:
            out.write("\t")
            cls = allCls[oId][marker]
            if cls is None:
                out.write("NA;")
            else:
                for r in range(len(cls)):
                    if r < len(com[oId]):
                        if cls[r] != com[oId][r]:
                            out.write("%s;" % cls[r])
                        elif not com[oId][r] is None:
                            out.write("--;")
                    else:
                        out.write("%s;" % cls[r])
        out.write("\n")
with open(snakemake.output.diffStat, "w") as sOut:
    for oId, data in stat.items():
        for mrk, mrkData in data.items():
            for oMrk, values in mrkData.items():
                sOut.write("%s\t%i\t%s\t%s\t%s\t%s\n" % (oId, size[oId], fungi[oId], mrk, oMrk, "\t".join(values)))
