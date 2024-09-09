from pylab import *
import read_sagecal as rt

cl = rt.getClusters(
    "/home/users/chege/theleap/leap/models/rescaled_new3c61.sky.txt.cluster",
    "/home/users/chege/theleap/leap/models/rescaled_new3c61.sky.txt",
)

output = "sky_bbs_DI.txt"
myf = open(output, "w")
myf.write(
    "format = Patch, Name, Type, Ra, Dec, I, SpectralIndex, ReferenceFrequency \n"
)
clsort = sorted(cl[0], key=lambda x: int(x["id"]))
for i in clsort:
    # if i['id']<=0:
    #        continue
    myf.write(
        "%s,,,%s,%s,,,\n"
        % (int(i["id"]), rt.radtohrstr(i["Ra"]), rt.radtodegstr(i["Dec"]))
    )
    for srcn in i["sources"].keys():
        src = i["sources"][srcn]
        myf.write(
            "%s,%s,POINT,%s,%s,%f,[%s],%f\n"
            % (
                "patch" + str(i["id"]),
                srcn,
                rt.radtohrstr(src["Ra"]),
                rt.radtodegstr(src["Dec"]),
                src["I"],
                ",".join(src["sp"]),
                src["freq0"],
            )
        )
        # myf.write("%s,%s,POINT,%s,%s,%f,[%s],%f\n"%(int(i['id']),srcn,rt.radtohrstr(src['Ra']),rt.radtodegstr(src['Dec']),src['I'],",".join(src['sp']),src['freq0']))
