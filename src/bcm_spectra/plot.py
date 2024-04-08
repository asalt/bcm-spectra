import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec


import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
import spectrum_utils.iplot as supi
from tqdm import tqdm
import altair as alt


def make_table_mpl(ax, table_data):
    ax.axis("tight")
    ax.axis("off")
    # table = ax.table(cellText=table_data.values, colLabels=table_data.columns, loc='center', cellLoc='center')
    table = ax.table(
        cellText=table_data.values,
        colLabels=table_data.columns,
        loc="center",
        cellLoc="center",
        fontsize=12,
        colWidths=[0.1, 0.1, 0.1],
    )
    return table


def make_table_altair(table_data):
    # not using this right now
    # Convert the DataFrame to a long format
    table_data_long = table_data.melt(var_name="Column", value_name="Text")

    # Creating a "table" using text marks
    text_chart = (
        alt.Chart(table_data_long)
        .mark_text(align="left", baseline="middle", dx=5)
        .encode(
            y=alt.Y("Column:N", axis=alt.Axis(title="")),
            x=alt.X("row_number:O", axis=alt.Axis(title=""), sort=None),
            text="Text:N",
            detail="Column:N",
        )
        .transform_window(row_number="row_number()")
        .transform_filter(
            alt.datum.Column
            != "Column"  # Optional filter if you have a 'Column' column
        )
    )
    return text_chart



def make_specutils_plot(msms, proforma_sequence='', additional_title=''):

    fig = plt.figure(figsize=(14, 6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    ax1 = fig.add_subplot(gs[0, 0])  # First row, first column
    ax2 = fig.add_subplot(gs[1, 0])  # Second row, first column
    sup.spectrum(msms, ax=ax1, grid=False)
    sup.mass_errors(msms, ax=ax2, plot_unknown=False)

    fig.suptitle(proforma_sequence, fontsize=16, weight="bold", color="black")
    fig.axes[0].title.set_text(additional_title)

    # table.set_fontsize(12)

    return fig