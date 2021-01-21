from smtplib import SMTP
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.header import Header
 
#邮件正文模板
SendHtml = """
<html>
    <head>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@4.5.0/dist/css/bootstrap.min.css" integrity="sha384-9aIt2nRpC12Uk9gS9baDl411NQApFmC26EwAOH8WgZl5MYYxFfc+NcPb1dKGj7Sk" crossorigin="anonymous">
    </head>
    <body background-color='red'>
    <p>Hi,<br>
        if you like this project?<br>
        <a href="http://www.realpython.com">NGSPipeDb</a> 
        Please star it on Github.
    </p>
    <p>{}</p>
    </body>
</html>
"""
def getAttachment(attachmentFilePath):  
    contentType, encoding = mimetypes.guess_type(attachmentFilePath)  
  
    if contentType is None or encoding is not None:  
        contentType = 'application/octet-stream'  
  
    mainType, subType = contentType.split('/', 1)  
    file = open(attachmentFilePath, 'rb')  
  
    if mainType == 'text':  
        attachment = MIMEText(file.read())  
    elif mainType == 'message':  
        attachment = email.message_from_file(file)  
    elif mainType == 'image':  
        attachment = MIMEImage(file.read(),_subType=subType)  
    elif mainType == 'audio':  
        attachment = MIMEAudio(file.read(),_subType=subType)  
    else:  
        attachment = MIMEBase(mainType, subType)  
    attachment.set_payload(file.read())  
    encode_base64(attachment)  
  
    file.close()

print ("输入发送邮箱服务器")
smtpserver = "smtp.sina.com"
print ("输入发送人邮箱用户/密码")
user = "ngspipedb@sina.com"
password = "8554ee2faf7409fe"
print ("输入发送人邮箱")
sender = "ngspipedb@sina.com"
print ("输入接收人邮箱")
receiver = ["wangjie191@mails.ucas.ac.cn"]
print ("输入接收人(抄送)邮箱") 
C_receiver = ["862730608@qq.com"]
print ("输入发送邮件主题")
subject = "Mail Python Test <验证>"
print ("编写Html类型邮件正文")
msg = MIMEMultipart("related")
msg["Subject"]=Header(subject,"utf-8")
msg['From'] = Header(user) #发件人
msg['to'] = Header(",".join(receiver))  #收件人（主送）
msg['Cc'] = Header(",".join(C_receiver))  #收件人（抄送）
msgText = MIMEText(SendHtml,"html","utf-8")
msg.attach(msgText) #添加邮件正文到邮件中
 
print ("输入发送附件")
attachmentFilePath = "/data/wangjie/mails/motif.png".encode("utf-8").decode("utf-8")
sendfile = open(attachmentFilePath,"rb").read()
addfile = MIMEText(sendfile,"base64","utf-8")
addfile["Content-Type"] = "application/octet-stream"
addfile["Content-Disposition"] = "attachment; filename={}".format(attachmentFilePath.split("\\")[-1])
msg.attach(addfile) #添加附件到邮件中
 
print ("连接发送邮件")
smtp = SMTP()
smtp.connect(smtpserver)
smtp.login(user, password)
smtp.sendmail(from_addr=sender, 
            to_addrs=receiver+C_receiver, 
            msg=msg.as_string())
smtp.quit()
