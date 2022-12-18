
import matplotlib.pyplot as plt
import demes
import demesdraw
import sys

#get the model id
model_id = sys.argv[1]

#load the model
graph = demes.load("{}.yaml".format(model_id))

#from tskit.dev/tutorials/viz.html
def size_max(graph):
    return max(max(epoch.start_size, epoch.end_size) for deme in graph.demes for epoch in deme.epochs)


w = 1.5 * size_max(graph)

positions = dict(ANC=0, B1=-w, B2=w, B1C1=-w-(w/2), B1C2=-w+(w/2), B2C1=w-(w/2),B2C2=w+(w/2),
                B1C1D1=-w-(1.5*w/2), B1C2D1=-w+(0.5*w/2), B2C1D1=w-(1.5*w/2),B2C2D1=w+(0.5*w/2),
                        B1C1D2=-w-(0.5*w/2), B1C2D2=-w+(1.5*w/2), B2C1D2=w-(0.5*w/2),B2C2D2=w+(1.5*w/2))


ax = demesdraw.tubes(graph, positions=positions, seed=42,log_time=True, labels="xticks-legend")
ax.set_title(model_id)

times=[epoch.end_time for deme in graph.demes for epoch in deme.epochs]
demenames=[ deme.name for deme in graph.demes ]


# get the x position of the text 
xpos = [positions[migration.source] + (positions[migration.dest] - positions[migration.source])/2 
                for migration in graph.migrations]

# get the y position of the text
ypos = [0.2*(graph.demes[demenames.index(migration.source)].epochs[-1].end_time + graph.demes[demenames.index(migration.dest)].epochs[0].start_time) 
                for migration in graph.migrations]

# get the labels of the text
labels = [f'{migration.rate:.2e}' for migration in graph.migrations]

for x,y,l in zip(xpos,ypos,labels):
        plt.text(x,y,l,fontsize=8)

plt.savefig("{}.png".format(model_id))


