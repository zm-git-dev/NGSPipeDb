#!/usr/bin/env python
# coding:utf-8
# @File : demo.py

import smtplib, sys
from email.mime.text import MIMEText


class Msmtp():
    def __init__(self, target, subject, content): # 收件人、标题、内容
        self.msg_from = '发件人QQ邮箱@qq.com'  # 邮件发送者
        self.password = '发件人QQ邮箱授权码'
        self.sender = smtplib.SMTP_SSL("smtp.qq.com", 465)
        self.msg_to = target.split(",")
        print self.msg_to
        self.subject = subject
        self.content = content

    def _login(self):
        self.sender.login(self.msg_from, self.password)

    def _msg(self):
        self.msg = MIMEText(self.content)  # 此处可选择文本格式或html等格式, 显示发送信息
        self.msg['Subject'] = self.subject
        self.msg['From'] = self.msg_from
        self.msg['To'] = ",".join(self.msg_to)

    def send_mail(self):
        try:
            self._login()
            self._msg()
            # sendmail 第二个参数，目的邮箱，参数类型 str 或者 list
            self.sender.sendmail(self.msg_from, self.msg_to, self.msg.as_string())
        except Exception, e:
            print u'邮件发送失败，原因：{}'.format( e)
        else:
            print u'邮件发送至 {} 成功！'.format(self.msg['To'])
        finally:
            self.sender.quit()

if __name__ == '__main__':
    # 收件人， 标题， 内容
    a = Msmtp(sys.argv[1], sys.argv[2], sys.argv[3])
    a.send_mail()
复制代码