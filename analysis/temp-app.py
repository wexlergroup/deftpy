import pandas as pd

eb_df = pd.read_csv("../data/features/Eb.csv")
vr_df = pd.read_csv("../data/features/Vr.csv")
bv_df = pd.read_csv("../data/features/Bv.csv")

from dash import dcc, html
from dash.dependencies import Input, Output, State
import dash

app = dash.Dash(__name__)

app.layout = html.Div([
    html.H2('Eb Dataset Lookup'),
    html.Div([
        html.Label('Search by Compound Name:'),
        dcc.Input(id='compound-name-input', type='text', placeholder='Enter compound name...'),
    ]),
    html.Div([
        html.Label('Search by Phase Info:'),
        dcc.Input(id='phase-info-input', type='text', placeholder='Enter phase info...'),
    ]),
    html.Div([
        html.Label('Minimum Value:'),
        dcc.Input(id='min-value-input', type='number', placeholder='Enter min value...'),
        html.Label('Maximum Value:'),
        dcc.Input(id='max-value-input', type='number', placeholder='Enter max value...'),
    ]),
    html.Button('Search', id='search-button', n_clicks=0),
    html.Div(id='query-output')
])

from dash.exceptions import PreventUpdate

@app.callback(
    Output('query-output', 'children'),
    [Input('search-button', 'n_clicks')],
    [State('compound-name-input', 'value'),
     State('phase-info-input', 'value'),
     State('min-value-input', 'value'),
     State('max-value-input', 'value')]
)
def update_output(n_clicks, compound_name, phase_info, min_value, max_value):
    if n_clicks < 1:
        raise PreventUpdate

    filtered_df = eb_df

    # Apply filters
    if compound_name:
        filtered_df = filtered_df[filtered_df['compound_name'].str.contains(compound_name, case=False, na=False)]
    if phase_info:
        filtered_df = filtered_df[filtered_df['phaseinfo'].str.contains(phase_info, case=False, na=False)]
    if min_value is not None:
        filtered_df = filtered_df[filtered_df['value'] >= min_value]
    if max_value is not None:
        filtered_df = filtered_df[filtered_df['value'] <= max_value]
    
    # Convert the filtered dataframe to a markdown table for display
    result_md = "```\n" + filtered_df.to_csv(index=False) + "\n```"
    
    return dcc.Markdown(result_md)

if __name__ == '__main__':
    app.run_server(debug=True)
