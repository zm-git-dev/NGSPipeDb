import sys
import smtplib, ssl
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.utils import formataddr

# ref: https://realpython.com/python-send-email/

# config
smtp_server = 'smtp.qq.com'
port = 465
sender_email = '296373256@qq.com'
receiver_email=sys.argv[1]
Authorization_code = 'rmcwksassrkgbjcf'


message = MIMEMultipart("alternative")
message["Subject"] = "inform from ngsPipe! work compeled!"
message["From"] = formataddr(["ngsPipe", sender_email])
message["To"] = formataddr(["", receiver_email])

# Create the plain-text and HTML version of your message

html = """\
<html>
    <body>
    <p>Hi,<br>
        if you like this project?<br>
        <a href="http://www.realpython.com">ngsPipe</a> 
        Please star it on Github.
    </p>
    </body>
</html>
"""

text = """\
rule1,rule2,stastics!"""

# Turn these into plain/html MIMEText objects
part1 = MIMEText(html, "html")
part2 = MIMEText(text, "plain",'utf-8')


# Add HTML/plain-text parts to MIMEMultipart message
# The email client will try to render the last part first
message.attach(part1)
message.attach(part2)


def mail():
    ret = True
    try:
        # Create secure connection with server and send email
        server = smtplib.SMTP_SSL(smtp_server, port)
        server.login(sender_email, Authorization_code)
        server.sendmail(
            sender_email, [receiver_email,], message.as_string())
        server.quit()
    except Exception:
        ret = False
    return ret
 
ret = mail()

if ret:
    print("邮件发送成功")
else:
    print("邮件发送失败")