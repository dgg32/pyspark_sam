import plotly.graph_objects as go
import json
import sys
inputfile = sys.argv[1]


data = json.load(open(inputfile))

fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = data["labels"],
      color = data["colors"],
      hovertemplate='Mutation value: %{value}<extra></extra>',
    ),
    link = dict(
      source = data["sources"], # indices correspond to labels, eg A1, A2, A1, B1, ...
      target = data["targets"],
      value =  data["values"]
  ))])

fig.update_layout(title_text="SAM mutations", font_size=10)
fig.show()