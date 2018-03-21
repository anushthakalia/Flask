from flask import Flask, render_template
# mapping multiple urls
app = Flask(__name__)

@app.route("/")
@app.route("/<user>")
def index(user=None):
	return render_template("multi_url.html",user= user)

if __name__=="__main__":
	app.run()