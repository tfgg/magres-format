def list_view(self):
    html = []

    html.append("<style>")
    html.append("table.efg, table.efg > tbody > tr > td, table.efg > tbody > tr { border:0; }")
    html.append("table.efg > tbody > tr { border-bottom:1px dashed #ccc; }")
    html.append("table.efg > tbody > tr:last-child { border-bottom:0; }")
    html.append("table.efg > tbody > tr > td > label { width:7em; color:#333; }")

    html.append("table.tensor { margin:0.5em; }")
    html.append("table.tensor, table.tensor > tbody > tr > td, table.tensor > tbody > tr { border:0; }")
    html.append("table.tensor > tbody > tr { border-bottom:1px solid #ccc; }")
    html.append("table.tensor > tbody > tr:last-child { border-bottom:0; }")
    html.append("table.tensor > tbody > tr > td { border-right:1px solid #ccc; padding:0.5em; }")
    html.append("table.tensor > tbody > tr > td:last-child { border-right:0; }")
    html.append("</style>")

    html.append("<table class='list'><tbody>")

    html.append(self[0]._repr_html_row_head_())
    for item in self:
      html.append(item._repr_html_row_())

    html.append("</tbody></table>")

    return "\n".join(html)

def tensor(V):
    html = ["<table class='tensor'>"]

    html.append("<tr><td>{:.3f}</td><td>{:.3f}</td><td>{:.3f}</td></tr>".format(*V[0,:]))
    html.append("<tr><td>{:.3f}</td><td>{:.3f}</td><td>{:.3f}</td></tr>".format(*V[1,:]))
    html.append("<tr><td>{:.3f}</td><td>{:.3f}</td><td>{:.3f}</td></tr>".format(*V[2,:]))

    html.append("</table>")

    return "\n".join(html)

def efg_row_head(self):
  return "<thead><tr><td>Atom</td><td>Q</td><td>Cq</td></tr></thead>"

def efg_row(self):
  return "<tr><td>{}</td><td>{:.3f}</td><td>{:.3f}</td></tr>".format(self.atom, self.atom.Q, self.Cq)

def efg(self):
    html = ["<h1>EFG on " + str(self.atom) + "</h1>"]

    html.append("<style>")
    html.append("table.efg, table.efg > tbody > tr > td, table.efg > tbody > tr { border:0; }")
    html.append("table.efg > tbody > tr { border-bottom:1px dashed #ccc; }")
    html.append("table.efg > tbody > tr:last-child { border-bottom:0; }")
    html.append("table.efg > tbody > tr > td > label { width:7em; color:#333; }")

    html.append("table.tensor { margin:0.5em; }")
    html.append("table.tensor, table.tensor > tbody > tr > td, table.tensor > tbody > tr { border:0; }")
    html.append("table.tensor > tbody > tr { border-bottom:1px solid #ccc; }")
    html.append("table.tensor > tbody > tr:last-child { border-bottom:0; }")
    html.append("table.tensor > tbody > tr > td { border-right:1px solid #ccc; padding:0.5em; }")
    html.append("table.tensor > tbody > tr > td:last-child { border-right:0; }")
    html.append("</style>")

    html.append("<table class='efg'><tbody>")

    html.append("<tr><td><label>Isotope</label></td><td>{}</td></tr>".format(self.atom.isotope))
    html.append("<tr><td><label>Q</label></td><td>{}</td></tr>".format(self.atom.Q))
    html.append("<tr><td><label>Cq</label></td><td>{:.3f} MHz</td></tr>".format(self.Cq))
    html.append("<tr><td><label>V</label></td><td>{}</td></tr>".format(tensor(self.V)))

    html.append("</tbody></table>")

    return "\n".join(html)

def ms_row_head(self):
  return "<thead><tr><td>Atom</td><td>ms iso</td></tr></thead>"

def ms_row(self):
  return "<tr><td>{}</td><td>{:.3f}</td></tr>".format(self.atom, self.iso)

def ms(self):
    html = ["<h1>MS on " + str(self.atom) + "</h1>"]

    html.append("<style>")
    html.append("table.ms, table.ms > tbody > tr > td, table.ms > tbody > tr { border:0; }")
    html.append("table.ms > tbody > tr { border-bottom:1px dashed #ccc; }")
    html.append("table.ms > tbody > tr:last-child { border-bottom:0; }")
    html.append("table.ms > tbody > tr > td > label { width:7em; color:#333; }")

    html.append("table.tensor { margin:0.5em; }")
    html.append("table.tensor, table.tensor > tbody > tr > td, table.tensor > tbody > tr { border:0; }")
    html.append("table.tensor > tbody > tr { border-bottom:1px solid #ccc; }")
    html.append("table.tensor > tbody > tr:last-child { border-bottom:0; }")
    html.append("table.tensor > tbody > tr > td { border-right:1px solid #ccc; padding:0.5em; }")
    html.append("table.tensor > tbody > tr > td:last-child { border-right:0; }")
    html.append("</style>")

    html.append("<table class='ms'><tbody>")

    html.append("<tr><td><label>Isotope</label></td><td>{}</td></tr>".format(self.atom.isotope))
    html.append("<tr><td><label>gamma</label></td><td>{:.3e}</td></tr>".format(self.atom.gamma))
    html.append("<tr><td><label>ms</label></td><td>{:.3f} ppm</td></tr>".format(self.iso))

    html.append("</tbody></table>")

    return "\n".join(html)

def isc(self):
    html = ["<h1>ISC between {:s} and {:s}</h1>".format(self.atom1, self.atom2)]

    html.append("<style>")
    html.append("table.isc, table.isc > tbody > tr > td, table.isc > tbody > tr { border:0; }")
    html.append("table.isc > tbody > tr { border-bottom:1px dashed #ccc; }")
    html.append("table.isc > tbody > tr:last-child { border-bottom:0; }")
    html.append("table.isc > tbody > tr > td > label { width:7em; color:#333; }")

    html.append("table.tensor { margin:0.5em; }")
    html.append("table.tensor, table.tensor > tbody > tr > td, table.tensor > tbody > tr { border:0; }")
    html.append("table.tensor > tbody > tr { border-bottom:1px solid #ccc; }")
    html.append("table.tensor > tbody > tr:last-child { border-bottom:0; }")
    html.append("table.tensor > tbody > tr > td { border-right:1px solid #ccc; padding:0.5em; }")
    html.append("table.tensor > tbody > tr > td:last-child { border-right:0; }")
    html.append("</style>")

    html.append("<table class='isc'><tbody>")

    html.append("<tr><td><label>Isotope 1</label></td><td>{}</td>".format(self.atom1.isotope))
    html.append("<td><label>Isotope 2</label></td><td>{}</td></tr>".format(self.atom2.isotope))
    html.append("<tr><td><label>Gamma 1</label></td><td>{}</td>".format(self.atom1.gamma))
    html.append("<td><label>Gamma 2</label></td><td>{}</td></tr>".format(self.atom2.gamma))
    html.append("<tr><td><label>K<sub>iso</sub></label></td><td>{:.3f}</td>".format(self.K_iso))
    html.append("<td><label>J<sub>iso</sub></label></td><td>{:.3f} Hz</td></tr>".format(self.J_iso))
    html.append("<tr><td><label>K</label></td><td colspan='3'>{}</td></tr>".format(tensor(self.K)))

    html.append("</tbody></table>")

    return "\n".join(html)

