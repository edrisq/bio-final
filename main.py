import datetime

from fasta.fasta import acc_to_ORFs
from fasta.align import needleman_wunsch, smith_waterman

from flask import Flask, render_template, request, jsonify

app = Flask(__name__)

@app.route('/')
def root():
    return render_template('index.html')

@app.route('/table_init', methods=['GET','POST'])
def table_init():
    seq1 = request.form['seq1']
    seq2 = request.form['seq2']

    return render_template('init.html', seq1=seq1, seq2=seq2)

@app.route('/table_fill', methods=['GET', 'POST'])
def table_fill():
    seq1 = request.form['seq1']
    seq2 = request.form['seq2']
    index = request.form['index']
    algo = request.form['algo']
    matrix = request.form['matrix']
    gap = request.form['gap']

    imax = int(index) // len(seq1) + 1
    jmax = int(index) % len(seq1) + 1

    if algo == 'nw':
        (f, sources) = needleman_wunsch(seq1, seq2, matrix, gap)
    else:
        (f, sources) = smith_waterman(seq1, seq2, matrix, gap)

    # return jsonify({'seq1' : seq1, 'seq2' : seq2, 'index' : index, 'algo' : algo, 'matrix' : matrix, 
    #                 'gap' : gap, 'f': f, 'sources': sources, 'imax': imax, 'jmax': jmax})
    return render_template('fill.html', seq1=seq1, seq2=seq2, index=index, algo=algo, 
                            matrix=matrix, gap=gap, f=f, sources=sources, imax=imax, jmax=jmax)

@app.route('/orf_finder')
def orf_root():
    return render_template('orf.html')

@app.route('/orf_results', methods=['GET','POST'])
def orf_finder():
    acc = request.form.get('acc') #if key doesn't exist, returns None
    ORFs = acc_to_ORFs(acc).split('\n')
    # print(ORFs)
    return render_template('results.html', acc=acc, ORFs=ORFs)

if __name__ == '__main__':
    # This is used when running locally only. When deploying to Google App
    # Engine, a webserver process such as Gunicorn will serve the app. This
    # can be configured by adding an `entrypoint` to app.yaml.
    # Flask's development server will automatically serve static files in
    # the "static" directory. See:
    # http://flask.pocoo.org/docs/1.0/quickstart/#static-files. Once deployed,
    # App Engine itself will serve those files as configured in app.yaml.
    app.run(host='127.0.0.1', port=8080, debug=True)