from flask import Flask, render_template

app= Flask(__name__)

@app.route('/')
@app.route('/<user>')
def index(user=None):
	return render_template("multi_url.html",user=user)

@app.route("/food")
def shopping():
	food=['cheese','burger','pizza']
	return render_template("food.html",food=food)

if __name__=="__main__":
	app.run()