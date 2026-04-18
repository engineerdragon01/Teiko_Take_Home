.PHONY: setup pipeline dashboard

setup:
	pip3 install -r requirements.txt

pipeline:
	python3 load_data.py
	python3 analysis.py

dashboard:
	python3 -m streamlit run app.py
