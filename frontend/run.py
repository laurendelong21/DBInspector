from os import path as pt
import time
import pandas as pd
from flask import Flask, request, redirect, render_template
from dbinspector.startup import CACHE, UNIPROT, REFSEQ
import dbinspector.parse
from dbinspector.compare import compare_entries, summary_statistics
from dbinspector.exceptions import QueryNotFoundError
from dbinspector.utils import determine_identifier_type, format_list_entry

UPLOAD_FOLDER = ''
ALLOWED_EXTENSIONS = {}

app = Flask(__name__)
app.secret_key = "someSecretKey;)" + time.strftime('%d%H')
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_PATH'] = 10 * 1024 * 1024  # Max 10MB


@app.route("/")
def home():
    return redirect('/summary')


@app.route("/summary")
def summary():
    """
    Route for the summary or home page. Checks if the necessary files exists, if not an button to generate them
    will be displayed. Otherwise the summary table will be shown.
    :return:
    """
    if pt.exists(CACHE) and pt.exists(pt.join(REFSEQ, 'refseq.json')) and pt.exists(pt.join(UNIPROT, 'uniprot.json')):
        try:
            summary: pd.DataFrame = summary_statistics()
            summary_html = summary.to_html(classes='data table table-striped', header="true", index=True, border=0,
                                           justify='left', na_rep=' ', table_id="results")
            return render_template('home.html', results=summary_html, parsed=True,
                                   current_time=time.strftime('%d.%m.%Y'))
        except Exception:
            return render_template('home.html', error="Summary could not be generated, try delete the cache and rerun!",
                                   current_time=time.strftime('%d.%m.%Y'))
    else:
        return render_template('home.html', message='Load databases!', current_time=time.strftime('%d.%m.%Y'))


@app.route("/comparison")
def comparison():
    """
    Route for the comparison page. With no input the lookup form is displayed.
    :return:
    """
    return render_template('comparison.html', current_time=time.strftime('%d.%m.%Y'))


@app.route('/populate', methods=['GET'])
def populate():
    """
    If necessary cache files do not exist, this REST route is called. It kicks of the parsing, mapping and
    saving of the databases.
    :return:
    """
    dbinspector.parse.parse_all()
    return redirect('/summary')


@app.route('/info', methods=['GET'])
def get_info():
    """
    The comparison form calls this REST route. Type of identifier is estimated and used to lookup and compare
    the entries in the cached databases.
    :return:
    """
    if request.method == 'GET':
        # print(request.args)
        symbol = request.args.get('identifier')
        query = determine_identifier_type(symbol)
        if query:
            try:
                query_res: pd.DataFrame = compare_entries(**query)
                query_res = query_res.applymap(lambda c: format_list_entry(c) if type(c) == list else c)
                res_html = query_res.to_html(classes='data table table-striped', header="true", index=True,
                                             border=0, justify='left', na_rep=' ', table_id="results")
                return render_template('comparison.html', results=res_html, query=symbol,
                                       current_time=time.strftime('%d.%m.%Y'))
            except QueryNotFoundError:
                return render_template('comparison.html', error='The queried identifier was not found!',
                                       current_time=time.strftime('%d.%m.%Y'))
        else:
            return render_template('comparison.html', error='Empty input!', current_time=time.strftime('%d.%m.%Y'))

    return render_template('comparison.html', current_time=time.strftime('%d.%m.%Y'))


if __name__ == '__main__':
    app.run(debug=True)
