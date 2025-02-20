#!/usr/bin/env python3
"""
Created on 2025-02-20 (Thu) 09:13:27

FlowKit Tutorial

Ref:
- https://flowkit.readthedocs.io/en/latest/notebooks/flowkit-tutorial-part01-sample-class.html


@author: I.Azuma
"""
# %%
import bokeh
from bokeh.plotting import show
bokeh.io.output_notebook()

import os
import flowkit as fk
from flowkit._utils import plot_utils

os.chdir('C:/github/omics-util/fcs_analysis')

# %% Part2
# サンプルインスタンスの作成とチャネルの確認
fcs_path = './data/Specimen_001_Ctrl_1_001.fcs'
sample = fk.Sample(fcs_path)
display(sample.channels)

# 分布の確認
p = sample.plot_histogram('FSC-H', source='raw', bins=256)
show(p)

# logicle transform
logicle_xform = fk.transforms.LogicleTransform(
    param_t=262144,
    param_w=0.5,
    param_m=4.5,
    param_a=0
)
sample.apply_transform(logicle_xform)

p = sample.plot_scatter(12, 10, source='xform', subsample=True)
show(p)

# asinh transform
asinh_xform = fk.transforms.AsinhTransform(
    param_t=262144,
    param_m=4.0,
    param_a=0.0
)
sample.apply_transform(asinh_xform)

p = sample.plot_scatter(12, 10, source='xform', subsample=True)
show(p)

# WSPBiexTransform (FlowJoのデフォルト値)
# Note: these are the default values in FlowJo (and in FlowKit, so you don't have to remember these specific values)
biex_xform = fk.transforms.WSPBiexTransform(
    max_value=262144.000029,
    positive=4.418540,
    width=-10,
    negative=0
)
sample.apply_transform(biex_xform)

p = sample.plot_scatter(12, 10, source='xform', subsample=True)
show(p)


# %% Part3
# xmlファイルで定義されたゲーティング情報を読み込む（ISACという基準に準拠しているらしい）
gml_path = './data/8_color_ICS.xml'
g_strat = fk.parse_gating_xml(gml_path)

text = g_strat.get_gate_hierarchy(output='ascii')
print(text)
display(g_strat.get_gate_ids())
display(g_strat.get_gate('Singlets'))

# %% Part4: 手動でゲーティング情報を追加する
import copy
import bokeh
from bokeh.plotting import show
import numpy as np

import flowkit as fk

bokeh.io.output_notebook()

# Load the sample
sample = fk.Sample("./data/Specimen_001_APAP_Low_1_003.fcs")
#display(sample.channels)

# Note: these are the default values in FlowJo (and in FlowKit, so you don't have to remember these specific values)

asinh_xform = fk.transforms.AsinhTransform(
    param_t=262144,
    param_m=1,
    param_a=0.0
)
sample.apply_transform(asinh_xform, include_scatter=True)  # 散乱光にも反映する場合はinclude_scatter=True
p = sample.plot_scatter('FSC-H', 'SSC-A', source='xform', subsample=True)
show(p)

chan_a_idx = sample.get_channel_index('FSC-H')
events_a = sample.get_channel_events(chan_a_idx, 'xform')
print(events_a.min(), events_a.max()) # >> (0.0407106872626841, 1.0000007612005928)


# Gating
g_strat = fk.GatingStrategy()  # Create a new GatingStrategy instance
g_strat.add_transform('asinh1', asinh_xform)

dim_a = fk.Dimension('FSC-H', range_max=1, transformation_ref="asinh1")
dim_b = fk.Dimension('SSC-A', range_max=1, transformation_ref="asinh1")

# polygon gate1
vertices = [
    (0.2, 2e-2),
    (0.5, 2e-2),
    (0.7, 0.2),
    (0.8, 0.6),
    (0.8, 1.2),
    (0.4, 1.2),
    (0.2, 0.3)
]

poly_gate = fk.gates.PolygonGate(
    'poly1',
    dimensions=[dim_a, dim_b],
    vertices=vertices
)
g_strat.add_gate(poly_gate, gate_path=('root',))
res = g_strat.gate_sample(sample)
res.report

# polygon gate2
vertices = [
    (0.6, 0.1),
    (0.6, 0.5),
    (0.9, 0.5),
    (0.9, 0.1)
]
poly_gate = fk.gates.PolygonGate(
    'poly2',
    dimensions=[dim_a, dim_b],
    vertices=vertices
)
g_strat.add_gate(poly_gate, gate_path=('root',))
res = g_strat.gate_sample(sample)
res.report

# boolean gate
gate_refs = [
    {
        'ref': 'poly1',
        'path': ('root',),
        'complement': False
    },
    {
        'ref': 'poly2',
        'path': ('root',),
        'complement': False
    }
]
bool_gate = fk.gates.BooleanGate('test', 'and', gate_refs)
res = g_strat.gate_sample(sample)
res.report

tmp_gate = res.get_gate_membership('poly1')
tmp_gate.sum()

# Finally, make a scatter plot with the gated events highlighted.
p = sample.plot_scatter('FSC-H', 'SSC-A', source='xform', highlight_mask=tmp_gate)
show(p)


# %%
#!/usr/bin/env python3
"""
Created on 2025-02-20 (Thu) 16:24:10

リンパ球のゲーティングの設定

@author: I.Azuma
"""
# %%
import bokeh
from bokeh.plotting import show
bokeh.io.output_notebook()

import os
import flowkit as fk
from flowkit._utils import plot_utils

os.chdir('C:/github/omics-util/fcs_analysis')

# %% 1. FSC-H vs SSC-A
sample = fk.Sample("./data/Specimen_001_APAP_Low_1_003.fcs")
display(sample.channels)

asinh_xform = fk.transforms.AsinhTransform(
    param_t=262144,
    param_m=1,
    param_a=0.0
)
sample.apply_transform(asinh_xform, include_scatter=True)
#p = sample.plot_scatter('FSC-H', 'SSC-A', source='xform', subsample=True)
#show(p)

g_strat = fk.GatingStrategy()  # Create a new GatingStrategy instance
g_strat.add_transform('asinh1', asinh_xform)

dim_a = fk.Dimension('FSC-H', range_min=0.0, range_max=2, transformation_ref="asinh1")
dim_b = fk.Dimension('SSC-A', range_min=0.0, range_max=2, transformation_ref="asinh1")

# polygon gating
vertices = [
    (0.2, 2e-2),
    (0.5, 2e-2),
    (0.7, 0.2),
    (0.9, 0.6),
    (0.9, 1.5),
    (0.4, 1.5),
    (0.2, 0.3)
]

poly_gate = fk.gates.PolygonGate(
    'fsc_ssc',
    dimensions=[dim_a, dim_b],
    vertices=vertices
)
g_strat.add_gate(poly_gate, gate_path=('root',))
res = g_strat.gate_sample(sample)
res.report

# %% 2. FSC-H vs FSC-A
#p = sample.plot_scatter('FSC-H', 'FSC-A', source='xform', subsample=True)
#show(p)

dim_a = fk.Dimension('FSC-H', range_max=1, transformation_ref="asinh1")
dim_b = fk.Dimension('FSC-A', range_max=1, transformation_ref="asinh1")

# polygon gating
vertices = [
    (0.05, 0),
    (0.1, 0),
    (0.9, 1.0),
    (0.9, 1.2),
    (0.05, 0.2),
]

poly_gate = fk.gates.PolygonGate(
    'fsc_fsc',
    dimensions=[dim_a, dim_b],
    vertices=vertices
)
g_strat.add_gate(poly_gate, gate_path=('root','fsc_ssc'))
res = g_strat.gate_sample(sample)
res.report

# %% 3. CD45 vs Dump
asinh_xform2 = fk.transforms.AsinhTransform(
    param_t=262144,
    param_m=4.0,
    param_a=0.0
)
sample.apply_transform(asinh_xform2, include_scatter=True)
try:
    g_strat.add_transform('asinh2', asinh_xform2)
except:
    pass

dim_a = fk.Dimension('PE-A', range_max=1, transformation_ref="asinh2")
dim_b = fk.Dimension('7-AAD-A', range_max=1, transformation_ref="asinh2")

# polygon gating
vertices = [
    (0.6, 0.5),
    (0.82, 0.5),
    (0.82, 0.75),
    (0.6, 0.75),
]

poly_gate = fk.gates.PolygonGate(
    'cd45_dump',
    dimensions=[dim_a, dim_b],
    vertices=vertices
)
g_strat.add_gate(poly_gate, gate_path=('root','fsc_ssc','fsc_fsc'))
res = g_strat.gate_sample(sample)
res.report

# %% CD3 vs NK
asinh_xform3 = fk.transforms.AsinhTransform(
    param_t=262144,
    param_m=3.3,
    param_a=0.0
)
sample.apply_transform(asinh_xform3, include_scatter=True)
p = sample.plot_scatter('FITC-A', 'APC-A', source='xform', subsample=True)
show(p)

g_strat.add_transform('asinh3', asinh_xform3)


quad_div1 = fk.QuadrantDivider(
    'cd3-div',
    'FITC-A',
    compensation_ref='uncompensated',
    transformation_ref="asinh3",
    values=[0.4]
)
quad_div2 = fk.QuadrantDivider(
    'nk-div',
    'APC-A',
    compensation_ref='uncompensated',
    transformation_ref="asinh3",
    values=[0.2]
)
quad_divs = [quad_div1, quad_div2]

# the 2 dividers above will be used to divide the space into 4 quadrants
quad_1 = fk.gates.Quadrant(
    quadrant_id='CD3pos-NKpos',
    divider_refs=['cd3-div', 'nk-div'],
    divider_ranges=[(0.4, None), (0.2, None)]
)
quad_2 = fk.gates.Quadrant(
    quadrant_id='CD3pos-NKneg',
    divider_refs=['cd3-div', 'nk-div'],
    divider_ranges=[(0.4, None), (None, 0.2)]
)
quad_3 = fk.gates.Quadrant(
    quadrant_id='CD3neg-NKpos',
    divider_refs=['cd3-div', 'nk-div'],
    divider_ranges=[(None, 0.4), (0.2, None)]
)
quad_4 = fk.gates.Quadrant(
    quadrant_id='CD3neg-NKneg',
    divider_refs=['cd3-div', 'nk-div'],
    divider_ranges=[(None, 0.4), (None, 0.2)]
)
quadrants = [quad_1, quad_2, quad_3, quad_4]

# We can now construct our QuadrantGate
quad_gate1 = fk.gates.QuadrantGate(
    'quadgate1',
    dividers=quad_divs,
    quadrants=quadrants
)

g_strat.add_gate(quad_gate1, gate_path=('root','fsc_ssc','fsc_fsc','cd45_dump'))
res = g_strat.gate_sample(sample)

# %% Visualization
# boolean gate
gate_refs = [
    {
        'ref': 'fsc_ssc',
        'path': ('root',),
        'complement': False
    },
    """
    {
        'ref': 'fsc_fsc',
        'path': ('root', 'fsc_ssc'),
        'complement': False
    },
    {
        'ref': 'cd45_dump',
        'path': ('root', 'fsc_ssc', 'fsc_fsc'),
        'complement': False
    },
    
    {
        'ref': 'CD3pos-NKpos',
        'path': ('root', 'fsc_ssc', 'fsc_fsc', 'cd45_dump'),
        'complement': False
    },
    """
]
#bool_gate = fk.gates.BooleanGate('tmp', 'and', gate_refs)
res = g_strat.gate_sample(sample)
res.report

tmp_gate = res.get_gate_membership('fsc_ssc')
asinh_xform = fk.transforms.AsinhTransform(
    param_t=262144,
    param_m=1,
    param_a=0.0
)
sample.apply_transform(asinh_xform, include_scatter=True)
p = sample.plot_scatter('FSC-H', 'SSC-A', source='xform', highlight_mask=tmp_gate)
show(p)

tmp_gate = res.get_gate_membership('fsc_fsc')
asinh_xform = fk.transforms.AsinhTransform(
    param_t=262144,
    param_m=1,
    param_a=0.0
)
sample.apply_transform(asinh_xform, include_scatter=True)
p = sample.plot_scatter('FSC-H', 'FSC-A', source='xform', highlight_mask=tmp_gate)
show(p)

tmp_gate = res.get_gate_membership('cd45_dump')
asinh_xform2 = fk.transforms.AsinhTransform(
    param_t=262144,
    param_m=4.0,
    param_a=0.0
)
sample.apply_transform(asinh_xform2, include_scatter=True)
p = sample.plot_scatter('PE-A', '7-AAD-A', source='xform', highlight_mask=tmp_gate)
show(p)

tmp_gate = res.get_gate_membership('CD3pos-NKpos')
asinh_xform3 = fk.transforms.AsinhTransform(
    param_t=262144,
    param_m=3.3,
    param_a=0.0
)
sample.apply_transform(asinh_xform3, include_scatter=True)
p = sample.plot_scatter('FITC-A', 'APC-A', source='xform', highlight_mask=tmp_gate)
show(p)

print(g_strat.get_gate_hierarchy(output='ascii'))

# %%
