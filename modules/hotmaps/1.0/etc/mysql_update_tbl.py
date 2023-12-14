#!/usr/bin/env python

# A simple script to update the character limits in the mupit_modbase MySQL mutations table

import argparse
import MySQLdb

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--host',
                        type=str, required=True,
                        help='MySQL host')
    parser.add_argument('-u', '--user',
                        type=str, required=True,
                        help='MySQL user')
    parser.add_argument('-p', '--passwd',
                        type=str, required=True,
                        help='MySQL password')
    parser.add_argument('-d', '--db',
                        type=str, required=True,
                        help='MySQL database')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)

def create_table(cursor, db):
    cursor.execute('drop table if exists mutations')
    cursor.execute('create table mutations (source varchar(10), tissue varchar(100), structure_id varchar(20), residues varchar(10), occurrence int) engine=innodb')
    db.commit()

def touch_output(outfile):
    complete = open(outfile, 'w')
    complete.close()

def main(args):

    db = MySQLdb.connect(host = args["host"],
                    user = args["user"],
                    passwd = args["passwd"],
                    db = args["db"])

    cursor = db.cursor()

    try:
        create_table(cursor, db)
        touch_output(args["output"])
    except MySQLdb.Error as e:
        db.rollback()
        print(f"Error occurred during table creation: {e}")
    finally:
        cursor.close()
        db.close()

if __name__ == '__main__':
    args = parse_args()
    main(args)