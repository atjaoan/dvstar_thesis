import plotly.express as px
import typer
from pathlib import Path
import pandas as pd 

app = typer.Typer()

cwd = Path(__file__).parent
csv_path = cwd / "tmp/benchmarks/test.csv"

@app.command()
def missses_vs_setsize():
    df = pd.DataFrame().read_csv(csv_path, dtype=str)
    #fig = px.line(df, x="year", y="lifeExp", title='Life expectancy in Canada')
    #fig.show()

if __name__ == "__main__":
    app()