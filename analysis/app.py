import pandas as pd
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
from dash.dash_table import DataTable


eb_df = pd.read_csv("../data/features/Eb.csv")
vr_df = pd.read_csv("../data/features/Vr.csv")
bv_df = pd.read_csv("../data/features/Bv.csv")

app = dash.Dash(__name__)

app.layout = html.Div([
    dcc.Tabs(id="tabs", value='tab-lookup', children=[
        dcc.Tab(label='Dataset Lookup', value='tab-lookup'),
        dcc.Tab(label='Crystal File Analysis', value='tab-crystal-analysis'),
        dcc.Tab(label='Data Trends', value='tab-trends'),
    ]),
    html.Div(id='tabs-content')
])

@app.callback(Output('tabs-content', 'children'),
              Input('tabs', 'value'))
def render_content(tab):
    if tab == 'tab-lookup':
        return html.Div([
            dcc.Dropdown(
                id='dataset-dropdown',
                options=[
                    {'label': 'Eb - Crystal Bond Dissociation Energy', 'value': 'eb'},
                    {'label': 'Vr - Crystal Reduction Potential', 'value': 'vr'},
                    {'label': 'Bv - Bond Valence Ratios', 'value': 'bv'},
                ],
                value='eb'
            ),
            html.Div(id='lookup-content')
        ])
    elif tab == 'tab-crystal-analysis':
        return html.Div([
            # Layout for uploading a crystal file and displaying analysis results
        ])
    elif tab == 'tab-trends':
        return html.Div([
            # Interface elements for exploring data trends
        ])

@app.callback(Output('lookup-content', 'children'),
              Input('dataset-dropdown', 'value'),
              prevent_initial_call=True)
def update_lookup(selected_dataset):
    if selected_dataset == 'eb':
        return html.Div([
            html.Label('Search by Formula:'),
            dcc.Input(id='eb-search', type='text', placeholder='Enter formula...'),
            html.Button('Search Eb', id='search-eb', n_clicks=0),
            html.Div(id='eb-results') 
        ])
    elif selected_dataset == 'vr':
        # Similar setup for 'vr', adjusting for dataset specifics
        return html.Div([
            html.Label('Search by Element:'),
            dcc.Input(id='vr-search', type='text', placeholder='Enter element...'),
            html.Button('Search Vr', id='search-vr', n_clicks=0),
            html.Div(id='vr-results')
        ])
    elif selected_dataset == 'bv':  # 'bv'
        # And for 'bv', also tailored to its specifics
        return html.Div([
            html.Label('Search by Site:'),
            dcc.Input(id='bv-search', type='text', placeholder='Enter site...'),
            html.Button('Search Bv', id='search-bv', n_clicks=0),
            html.Div(id='bv-results')
        ])
    
@app.callback(
    Output('eb-results', 'children'),
    Input('search-eb', 'n_clicks'),
    State('eb-search', 'value'),
    prevent_initial_call=True 
)
def search_eb(n_clicks, search_value):
    if n_clicks > 0:
        
        filtered_df = eb_df[eb_df['formula'].str.contains(search_value, case=False, na=False)]
    
        return DataTable(
            id='eb-table',
            columns=[{"name": i, "id": i} for i in filtered_df.columns],
            data=filtered_df.to_dict('records'),
            style_table={'overflowX': 'auto'}, 
            filter_action="native",  
            sort_action="native", 
            sort_mode="multi",
            page_action="native",
            page_size=10,
        )
    else:
        return ''
    
if __name__ == '__main__':
    app.run_server(debug=True)
