from flask import Flask, render_template

app= Flask(__name__)

@app.route('/profile/<name>')  #not a good practise to pass html through return , html files usualyy are
#in the templates folder and the CSS files in the  static folder
def profile(name):
	return render_template("profile.html",name=name)

if __name__=="__main__":
	app.run(debug=True)