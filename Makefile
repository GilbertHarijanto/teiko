PYTHON ?= python3
VENV := .venv
PY := $(VENV)/bin/python

.PHONY: setup pipeline dashboard

setup:
	$(PYTHON) -m venv $(VENV)
	$(PY) -m pip install -r requirements.txt

pipeline:
	$(PY) load_data.py
	$(PY) frequency.py
	MPLBACKEND=Agg $(PY) response_analysis.py
	$(PY) subset_analysis.py

dashboard:
	-@lsof -ti tcp:8501 | xargs kill 2>/dev/null || true
	$(PY) -m streamlit run dashboard.py --server.headless true --server.address 0.0.0.0 --server.port 8501
