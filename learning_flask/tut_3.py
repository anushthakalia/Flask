from flask import Flask, request

app= Flask(__name__)

@app.route('/')
def index():
	return 'Method used: %s' % request.method

#when we request a URL, its using GET.It is the default request method
#But when we are submitting a form or something we use POST
@app.route('/user',methods=['GET','POST'])# the possible request methods this page can execute
def user():
	if request.method=='POST':
		return "its POST"
	else:
		return "its GET"

if __name__=="__main__":
	app.run(debug=True)


#try the post man app for testing out POST and GET