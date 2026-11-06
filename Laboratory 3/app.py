import sys

import pandas as pd
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import (
    QApplication, QComboBox, QGridLayout, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow,
    QPushButton, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

import analysis


class Worker(QThread):
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, function, *args):
        super().__init__()
        self.function, self.args = function, args

    def run(self):
        try:
            self.done.emit(self.function(*self.args))
        except Exception as error:
            self.failed.emit(str(error))


class Chart(FigureCanvasQTAgg):
    def __init__(self):
        super().__init__(Figure(figsize=(6, 4), layout="tight"))

    def bars(self, series, title, horizontal=False):
        self.figure.clear()
        axes = self.figure.add_subplot()
        if horizontal:
            axes.barh(series.index[::-1], series.values[::-1], color="#4c78a8")
        else:
            axes.bar(series.index, series.values, color="#4c78a8")
            axes.tick_params(axis="x", labelrotation=90, labelsize=7)
        axes.set_title(title)
        self.draw()


def make_table():
    table = QTableWidget()
    table.verticalHeader().hide()
    table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    return table


def fill_table(table, frame):
    table.setRowCount(len(frame))
    table.setColumnCount(len(frame.columns))
    table.setHorizontalHeaderLabels(list(frame.columns))
    for row, values in enumerate(frame.itertuples(index=False)):
        for column, value in enumerate(values):
            table.setItem(row, column, QTableWidgetItem(str(value)))
    header = table.horizontalHeader()
    header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
    table.setWordWrap(True)
    table.resizeRowsToContents()


class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Virus Diet Explorer")
        self.resize(1200, 800)
        self.matches = pd.DataFrame()
        self.worker = None

        self.query = QLineEdit(placeholderText="Virus name, e.g. influenza, coronavirus, Wanjuvirus PP2")
        self.search_button = QPushButton("Search")
        self.genomes = QComboBox()
        self.analyze_button = QPushButton("Analyze", enabled=False)
        self.status = QLabel("Type a virus name and press Search. Lab proof of concept, not medical advice.")
        self.metrics = QLabel()

        self.top_codons = make_table()
        self.codon_chart = Chart()
        self.amino_chart = Chart()
        self.summary = QLabel(wordWrap=True)
        self.eat = make_table()
        self.avoid = make_table()

        codon_tab = QWidget()
        codon_layout = QHBoxLayout(codon_tab)
        codon_layout.addWidget(self.top_codons, 1)
        codon_layout.addWidget(self.codon_chart, 2)

        diet_tab = QWidget()
        diet_layout = QGridLayout(diet_tab)
        diet_layout.addWidget(self.summary, 0, 0, 1, 2)
        diet_layout.addWidget(QLabel("<h3>✅ Eat</h3>"), 1, 0)
        diet_layout.addWidget(QLabel("<h3>⛔ Avoid</h3>"), 1, 1)
        diet_layout.addWidget(self.eat, 2, 0)
        diet_layout.addWidget(self.avoid, 2, 1)

        self.tabs = QTabWidget()
        self.tabs.addTab(diet_tab, "Diet")
        self.tabs.addTab(codon_tab, "Codons")
        self.tabs.addTab(self.amino_chart, "Amino acids")

        search_row = QHBoxLayout()
        search_row.addWidget(self.query, 1)
        search_row.addWidget(self.search_button)
        genome_row = QHBoxLayout()
        genome_row.addWidget(self.genomes, 1)
        genome_row.addWidget(self.analyze_button)

        root = QWidget()
        layout = QVBoxLayout(root)
        layout.addLayout(search_row)
        layout.addLayout(genome_row)
        layout.addWidget(self.status)
        layout.addWidget(self.metrics)
        layout.addWidget(self.tabs, 1)
        self.setCentralWidget(root)

        self.query.returnPressed.connect(self.search)
        self.search_button.clicked.connect(self.search)
        self.analyze_button.clicked.connect(self.analyze)

    def run(self, message, function, *args, then):
        self.status.setText(message)
        self.search_button.setEnabled(False)
        self.analyze_button.setEnabled(False)
        self.worker = Worker(function, *args)
        self.worker.done.connect(then)
        self.worker.failed.connect(lambda error: self.status.setText(f"Error: {error}"))
        self.worker.finished.connect(lambda: self.search_button.setEnabled(True))
        self.worker.finished.connect(lambda: self.analyze_button.setEnabled(not self.matches.empty))
        self.worker.start()

    def search(self):
        if self.query.text().strip():
            self.run("Searching NCBI...", analysis.search, self.query.text(), then=self.show_matches)

    def show_matches(self, matches):
        self.matches = matches
        self.genomes.clear()
        for row in matches.itertuples():
            self.genomes.addItem(" · ".join(part for part in (row.organism, row.strain, row.accession) if part))
        self.status.setText(f"{len(matches)} genomes found." if len(matches) else "No genomes found on NCBI for that name.")

    def analyze(self):
        genome = self.matches.iloc[self.genomes.currentIndex()]
        message = "Downloading the genome and the USDA food database (first run only), then analyzing..."
        self.run(message, analysis.analyze, genome["accession"], genome["ftp"], then=self.show_results)

    def show_results(self, result):
        codons = result["codons"]
        total = sum(codons.values())
        self.status.setText("Done.")
        self.metrics.setText(
            f"<b>Genome length:</b> {result['length']:,} nt &nbsp; <b>Segments:</b> {result['segments']} "
            f"&nbsp; <b>Codons counted:</b> {total:,}"
        )

        top = pd.DataFrame(codons.most_common(10), columns=["Codon", "Count"])
        top["Frequency"] = (top["Count"] / total * 100).round(2).astype(str) + "%"
        top["Amino acid"] = top["Codon"].map(lambda c: analysis.AMINO_ACID_NAMES[analysis.AMINO_ACIDS[c]])
        fill_table(self.top_codons, top)
        self.codon_chart.bars(pd.Series(codons).sort_values(ascending=False), "All 64 codons")
        self.amino_chart.bars(pd.Series(result["amino_acids"]).sort_values(ascending=False), "Amino acid profile", horizontal=True)

        self.summary.setText(
            f"This genome's codons favour <b>{', '.join(result['favourite'])}</b>. Foods from the USDA database "
            "whose protein has the most similar amino-acid mix are listed under Avoid, and the least similar under "
            f"Eat. Only foods with at least {analysis.MIN_PROTEIN} g of amino acids per 100 g are compared."
        )
        fill_table(self.eat, result["eat"])
        fill_table(self.avoid, result["avoid"])


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec())
