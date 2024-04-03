import base64
import re
import dash
from dash import html, dcc, Input, Output, State, callback
from io import StringIO
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import Poscar
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
from crystal_analysis import Crystal

# initialize Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div([
    dbc.Container([
        dbc.Row([
            dbc.Col(html.H1("Crystal Feature Dashboard"), className="mb-2")
        ]),
        dbc.Row([
            dbc.Col(html.H6("Upload a crystal file to analyze"), className="mb-4")
        ]),
        dbc.Row([
            dbc.Col(dcc.Upload(
                id='upload-data',
                children=html.Div(['Drag and Drop or ', html.A('Select a File')]),
                className="upload-div"  
            ), className="mb-4"),
        ]),
        dbc.Row([
            dbc.Col(dcc.Graph(id='eb-graph'), className="graph mb-4"), 
            dbc.Col(dcc.Graph(id='vr-graph'), className="graph mb-4"),  
            dbc.Col(dcc.Graph(id='bv-graph'), className="graph mb-4"), 
        ])
    ], className="container")  
])

def parse_contents(contents, filename):
    content_type, content_string = contents.split(',')
    
    decoded = base64.b64decode(content_string)
    
    if filename.lower().endswith('.cif'):
        content_string = decoded.decode('utf-8')
        parser = CifParser(StringIO(content_string))
        structure = parser.parse_structures(primitive=True)[0]
        crystal = Crystal(pymatgen_structure=structure)
    elif filename.lower().endswith('.vasp') or filename.lower().endswith('.poscar'):
        content_string = decoded.decode('utf-8')
        poscar = Poscar.from_string(content_string)
        crystal = Crystal(pymatgen_structure=poscar.structure)
    else:
        raise ValueError("File format not supported. Please upload a CIF or POSCAR file.")
    
    eb = crystal.bond_dissociation_enthalpies
    vr = crystal.reduction_potentials
    bv = crystal.bond_valence_sum_ratios
    
    return {
        'Eb': eb, 
        'Vr': vr, 
        'Bv': bv,
    }

@callback(
    [Output('eb-graph', 'figure'),
     Output('vr-graph', 'figure'),
     Output('bv-graph', 'figure')],
    [Input('upload-data', 'contents')],
    [State('upload-data', 'filename')]
)
def update_graph(contents, filename):
    if contents is not None:
        data = parse_contents(contents, filename)
        eb_trace = go.Scatter(x=list(range(len(data['Eb']))), y=data['Eb'], mode='lines+markers')
        vr_trace = go.Scatter(x=list(range(len(data['Vr']))), y=data['Vr'], mode='lines+markers')
        bv_trace = go.Scatter(x=list(range(len(data['Bv']))), y=data['Bv'], mode='lines+markers')

        eb_graph = {'data': [eb_trace], 'layout': go.Layout(title='Eb Trend')}
        vr_graph = {'data': [vr_trace], 'layout': go.Layout(title='Vr Trend')}
        bv_graph = {'data': [bv_trace], 'layout': go.Layout(title='Bv Trend')}

        return eb_graph, vr_graph, bv_graph
    else:
        # If no data has been uploaded return empty graphs
        return {}, {}, {}

if __name__ == '__main__':
    app.run_server(debug=True)
