from flask import Flask

app= Flask(__name__)

@app.route('/')
def index():
	return 'This is the homepage'

@app.route('/about')
def about():
	return "<h2>We are a very old company</h2>"

@app.route('/profile/<username>') # strings can be passed simply
def profile(username):
	return "<h2>Hey there %s</h2>" % username

@app.route('/post/<int:post_id>') #this is the way you pass integers
def post(post_id):
	return "<h2>Post ID is %s</h2>" % post_id

if __name__=="__main__":
	app.run(debug=True)